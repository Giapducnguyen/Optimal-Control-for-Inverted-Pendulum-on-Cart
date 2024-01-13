% Draft

close all; clear all; clc

%% Model data
M = 2;                  % [kg]
m = 0.1;                % [kg]
l = 0.5;                % [m]
g = 9.81;               % [m/s^2]
Ts = 0.01;              % [s]
Tspan = 0:Ts:30;        % [s]
Nsim = length(Tspan); 
x0 = [0;0;0;0];         % [rad; rad/s; m; m/s]
xr = 1;                 % [m]

u_max = 30; u_min = -u_max;                         % Control Constraints
x3_max = 1.05; x3_min = -0.55;                      % Cart Position Constraints
theta_max = deg2rad(15); theta_min = -theta_max;    % Rod Angle Constraints

%% Design Parameters
Np = 1/(Ts);                % Prediction Horizon
Qe = diag([1, 1]);          % Tracking Penalty
R = 0.5;                    % Control Increment Penalty

modeling_option = 0;        % Nonlinear Dynamics Modeling Option: whether to use MATLAB symbolic toolbox
terminalCost_option = 1;    % MPC Option: whether to use terminal cost
discretization_option = 0;  % Discretization Option: whether to use c2d or derive manually
ode45_option = 0;           % Simulation Option: whether to use ode45 or manual RK4

%% Nonlinear Dynamic Modeling
if modeling_option == 1 % derivation by solving symbolic expression (MATLAB)
    syms M m l g ddx theta dtheta ddtheta ut dt

    f_sym = [(M+m)*ddx - m*l*sin(theta)*dtheta^2 + m*l*cos(theta)*ddtheta - ut - dt;
            m*ddx*cos(theta) + m*l*ddtheta - m*g*sin(theta) - dt*cos(theta)];
    
    sol = solve(f_sym,[ddtheta, ddx]);
    ddtheta = simplify(sol.ddtheta);
    ddx = simplify(sol.ddx);

    nonlDyns = @(t,x,u,d) ([x(2); ...
              double(subs(ddtheta,[theta; dtheta; ut; dt],[x(1); x(2); u; d])); ...
              x(4); ...
              double(subs(ddx,[theta; dtheta; ut; dt],[x(1); x(2); u; d]))]);
else % manual derivation
    nonlDyns = @(t,x,u,d) ([x(2); ...
              (m^2*g*sin(x(1)) - m*cos(x(1))*u + M*cos(x(1))*d + M*m*g*sin(x(1))...
              - m^2*l*x(2)^2*sin(x(1))*cos(x(1)))/(m*l*(-m*cos(x(1))^2+M+m)); ...
              x(4); ...
              (m*l*sin(x(1))*x(2)^2 - d*cos(x(1))^2 - m*g*sin(x(1))*cos(x(1))...
              + d + u)/(-m*cos(x(1))^2+M+m)]);
end

%% Nonlinear Dynamics Propagation by Runge-Kutta 4th order
k1 = @(t,x,u,d) (nonlDyns(t,x,u,d));
k2 = @(t,x,u,d) (nonlDyns(t,x + k1(t,x,u,d)*Ts/2,u,d));
k3 = @(t,x,u,d) (nonlDyns(t,x + k2(t,x,u,d)*Ts/2,u,d));
k4 = @(t,x,u,d) (nonlDyns(t,x + k3(t,x,u,d)*Ts,u,d));
f_rk4 = @(t,x,u,d) (x + (Ts/6)*(k1(t,x,u,d) + 2*k2(t,x,u,d) + 2*k3(t,x,u,d) + k4(t,x,u,d)));

%% Linearized Dynamics
% x = [x1; x2; x3; x4] = [theta; dtheta; x; dx]
Ac = [0,              1, 0, 0;
      (M+m)*g/(M*l),  0, 0, 0;
      0,              0, 0, 1;
      -m*g/M,         0, 0, 0];
Buc = [0; -1/(M*l); 0; 1/M];
Bdc = [0;  1/(m*l); 0; 0];
Cc =  [1, 0, 0, 0;
       0, 0, 1, 0];
Dc =  0;

[nx,nu] = size(Buc);
ny = size(Cc,1);
nr = ny;

%% Discretization of linearized Dynamics
if discretization_option == 1 % discretization using MATLAB
    sysc = ss(Ac, [Buc, Bdc], Cc, Dc);
    d_method = 'zoh'; % 'tustin' or 'zoh'
    sysd = c2d(sysc,Ts,d_method);
    Ad = sysd.A;
    Bd_ext = sysd.B; Bud = Bd_ext(:,1:nu); Bdd = Bd_ext(:,nu+1:end);
    Cd = sysd.C;
    Dd = sysd.D;
else % manual derivation
    N_iter = 100;
    [Ad,Bud,Bdd] = discretization(Ac,Buc,Bdc,Ts,N_iter);
    Cd = Cc; Dd = Dc;
end

%% Extended dynamics
% % extended-state: x_ext = [dx; e; x-; u-]
A_ext = [Ad,            zeros(nx,nr),   zeros(nx),      zeros(nx,nu);
         Cd*Ad,         eye(nr),        zeros(nr,nx),   zeros(nr,nu);
         eye(nx),       zeros(nx,nr),   eye(nx),        zeros(nx,nu);
         zeros(nu,nx),  zeros(nu,nr),   zeros(nu,nx),   eye(nu)];
B_ext = [Bud; Cd*Bud;   zeros(nx,nu);   eye(nu)];

%% Control Constraints
% Polyhedral representation of input constraints: Au * u <= bu
Au = [1;        -1];
bu = [u_max;    -u_min];

Au_ext = Au*[zeros(nu,nx), zeros(nu,nr), zeros(nu,nx), eye(nu)];
bu_ext = bu;

%% State Constraints
% Polyhedral representation of state constraints: Ax * x <= bx
Ax = [ 1, 0, 0, 0;
      -1, 0, 0, 0;
       0, 0, 1, 0;
       0, 0,-1, 0];
bx = [theta_max; -theta_min; x3_max; -x3_min];

Ax_ext = Ax*[zeros(nx), zeros(nx,nr), eye(nx), zeros(nx,nu)];
bx_ext = bx;


%% MPC Data
Q_ext = blkdiag(zeros(nx), Qe, zeros(nx), zeros(nu));
R_ext = R;

%% Terminal Penalty Computation
Af = [Ad,        zeros(nx,nr);
      Cd*Ad,     eye(nr)];
Bf = [Bud; Cd*Bud];
Qf = blkdiag(zeros(nx), Qe);
Rf = R;

if terminalCost_option == 1
    [~,P1,~] = dlqr(Af,Bf,Qf,Rf,[]);
    % [P1,~,~] = idare(Af,Bf,Qf,Rf,[],[]);
    P = blkdiag(P1, zeros(nx), zeros(nu));
else
    P = Q_ext;
end

%% Stacked constraint over prediction horizon
Uad = admissibleInputs(A_ext,B_ext,Np,Au_ext,bu_ext,Ax_ext,bx_ext);

%% Quadratic programming matrices
Q_bar = blkdiag(kron(eye(Np-1),Q_ext),P);
R_bar = kron(eye(Np),R_ext);
[A_bar,B_bar] = liftedDynamics(A_ext,B_ext,Np);
H = B_bar'*Q_bar*B_bar + R_bar;    % H = (H+H')/2;
options = mpcActiveSetOptions;
iA0 = false(size(Uad.b));

%% Main simulation loop
Xlog = nan(nx, Nsim+1);
Xlog(:,1) = x0;
Ucomplog = nan(nu, Nsim);
Uactlog = nan(nu, Nsim);
Elog = zeros(nr, Nsim+1);
infeasibleOCP = 0;
dUcomplog = nan(nu,Nsim);
Xrlog = nan(1,Nsim);

for k = 1:Nsim
    % % Reference Update
    if k > (10/Ts)
        xr = -0.5;
    end
    r = [0; xr];
    Xrlog(:,k) = xr; 
    % % Disturbance Update
    % a nudge (impulse disturbance):  (k >= (20/Ts)) && (k <= (20/Ts)+5)
    % constant force:   k > (20/Ts)
    if (k >= (20/Ts)) && (k <= (20/Ts)+10)
        dk = 1;
    else
        dk = 0;
    end
    % % Current output
    yk = Cd*Xlog(:,k);

    % % State incremental computation
    if k == 1
        delta_xk = zeros(nx,1);
    else
        delta_xk = Xlog(:,k) - Xlog(:,k-1);
    end

    % % Current error computation
    ek = yk - r;
    Elog(:,k) = ek;

    % % Previous state
    if k == 1
        xk_1 = x0;
    else
        xk_1 = Xlog(:,k-1);
    end

    % % Previous control
    if k == 1
        uk_1 = 0;
    else
        uk_1 = Uactlog(:,k-1);
    end

    % % Extended state
    x_ext = [delta_xk; ek; xk_1; uk_1];
    
    % % Control computation
    % linear cost of QP computation
    q = B_bar'*Q_bar*A_bar*x_ext;

    % Control increment
    [dU,exitflag,iA0,~] = mpcActiveSetSolver(H,q,Uad.A,Uad.b-Uad.B*x_ext,Uad.Ae,Uad.be,iA0,options);
    % [dU,~,exitflag,~] = quadprog(H,q,Uad.A,Uad.b-Uad.B*x_ext);
    if exitflag <= 0
        delta_uk = 0; 
        infeasibleOCP = infeasibleOCP + 1;
    else
        delta_uk = dU(1,nu);
    end
    dUcomplog(:,k) = delta_uk;
    % (computed) current control
    uk = uk_1 + delta_uk;
    Ucomplog(:,k) = uk;
    % actual (saturated) current control
    Uactlog(:,k) = max(u_min,min(uk,u_max));

    % % Update system
    if ode45_option == 1
        [t,xk1] = ode45(@(t,x) nonlDyns(0,x,Uactlog(:,k),dk), [(k-1)*Ts, (k)*Ts], Xlog(:,k));
        Xlog(:,k+1) = xk1(end,:)';
    else
        Xlog(:,k+1) = f_rk4(0, Xlog(:,k), Uactlog(:,k), dk);
    end
end

%% Plots
lw = 1;
figure(1)
fig1 = tiledlayout(2,1);

nexttile
hold on;
plot(Tspan, Xlog(1,1:end-1),'LineStyle','-','LineWidth',lw,'Color','b');
plot(Tspan, theta_min*ones(size(Tspan)),'LineStyle','--','LineWidth',lw,'Color','k');
plot(Tspan, theta_max*ones(size(Tspan)),'LineStyle','--','LineWidth',lw,'Color','k');
legend({'LQMPC','Constraints'},'Interpreter','latex');
ylabel('$\theta (\mathrm{rad})$','Interpreter','latex');
title('(a) Rod Angle','Interpreter','latex')
ylim([-0.3 0.3]); yticks([-0.2618 0 0.2618]);
box on; grid on;
ax = gca; ax.FontSize = 14;

nexttile
hold on; box on; grid on;
plot(Tspan, Xlog(3,1:end-1),'LineStyle','-','LineWidth',lw,'Color','b');
plot(Tspan,Xrlog,'LineStyle','--','LineWidth',lw,'Color','k');
plot(Tspan, x3_min*ones(size(Tspan)),'LineStyle','--','LineWidth',0.5,'Color','r');
plot(Tspan, x3_max*ones(size(Tspan)),'LineStyle','--','LineWidth',0.5,'Color','r');
legend({'LQMPC','Reference'},'Interpreter','latex');
ylabel('$x (\mathrm{m})$','Interpreter','latex');
xlabel('Time (s)','Interpreter','latex');
title('(b) Cart Position','Interpreter','latex')
ylim([-0.7 1.2]); yticks([-0.5 0 1]);
ax = gca; ax.FontSize = 14;

fig1.TileSpacing = 'compact';
fig1.Padding = 'compact';

% % % % ------------------------------------------------- % % % %

figure(2)
fig2 = tiledlayout(2,1);

nexttile
stairs(Tspan,dUcomplog,'LineStyle','-','LineWidth',lw,'Color','b');
ylabel('$\delta u (\mathrm{N})$','Interpreter','latex');
title('(a) Incremental Control Force','Interpreter','latex');
ylim([-5 20]); yticks([-3 0 18]);
box on; grid on;
ax = gca; ax.FontSize = 14;

nexttile
hold on; box on; grid on;
stairs(Tspan,Ucomplog,'LineStyle','-','LineWidth',lw,'Color','b');
stairs(Tspan,Uactlog,'LineStyle','--','LineWidth',lw,'Color','r');
stairs(Tspan,u_max*ones(size(Tspan)),'LineStyle','--','LineWidth',lw,'Color','k');
stairs(Tspan,u_min*ones(size(Tspan)),'LineStyle','--','LineWidth',lw,'Color','k');
ylabel('$u (\mathrm{N})$','Interpreter','latex');
xlabel('Time (s)','Interpreter','latex');
legend({'Computed','Saturated','Constraints'},'Interpreter','latex');
title('(b) Control Force','Interpreter','latex')
ylim([-34 34]); yticks([-30 0 30]);
ax = gca; ax.FontSize = 14;

fig2.TileSpacing = 'compact';
fig2.Padding = 'compact';

%% Function helper 1: Discretization - zoh
function [Ad,Bud,Bdd] = discretization(Ac,Buc,Bdc,Ts,N_iter)
if det(Ac) ~= 0
    Ad = expm(Ac*ts);
    Bud = Ac\(Ad-eye(size(Ad)))*Buc;
    Bdd = Ac\(Ad-eye(size(Ad)))*Bdc;
else
    Ad = zeros(size(Ac));
    Bud = zeros(size(Buc));
    Bdd = zeros(size(Bdc));
    
    for k=0:N_iter
        Ad = Ad + (1/factorial(k))*(Ac*Ts)^k;
        Bud = Bud + (1/factorial(k+1))*Ac^(k)*Ts^(k+1)*Buc;
        Bdd = Bdd + (1/factorial(k+1))*Ac^(k)*Ts^(k+1)*Bdc;
    end
end
end

%% Function helper 2: constraint transformation to quadratic programming

function Uad = admissibleInputs(A,B,Np,Gu,gu,Fx,fx)

[A_bar,B_bar] = liftedDynamics(A,B,Np);

Uad.A = [kron(eye(Np),Fx)*B_bar;
         kron(eye(Np),Gu)*B_bar];
Uad.b = [kron(ones(Np,1),fx);
         kron(ones(Np,1),gu)];
Uad.B = [kron(eye(Np),Fx)*A_bar;
         kron(eye(Np),Gu)*A_bar];
Uad.Ae = zeros(0,size(B_bar,2));
Uad.be = zeros(0,1);
end

%% Function helper 3: Stacking dynamics over prediction horizon
function [A_bar,B_bar] = liftedDynamics(A,B,Np)

% A_bar = cell2mat(cellfun(@(x)A^x,num2cell((1:Np)'),'UniformOutput',false));
% B_bar = tril(cell2mat(cellfun(@(x)A^x,num2cell(toeplitz(0:Np-1)),'UniformOutput',false)))*kron(eye(Np),B);

% % or using cell

A_bar = cell(Np, 1);
B_bar = cell(Np,Np);
b0 = zeros(size(B));
for i = 1:Np
    A_bar{i} = A^i;
    for j = 1:Np
        if i >= j
            B_bar{i,j} = A^(i-j)*B;
        else
            B_bar{i,j} = b0;
        end
    end
end
A_bar = cell2mat(A_bar);
B_bar = cell2mat(B_bar);
end