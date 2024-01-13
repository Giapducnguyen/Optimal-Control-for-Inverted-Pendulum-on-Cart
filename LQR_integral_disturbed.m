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
u_max = 30; u_min = -u_max;
theta_max = deg2rad(15); theta_min = -theta_max;

%% Design Parameters
Q = diag([10, 1, 100, 1, 1]);   % Tracking Penalty
R = 2;                          % Control Penalty

modeling_option = 0;        % Nonlinear Dynamics Modeling Option: whether to use MATLAB symbolic toolbox (1) or manual derivation
discretization_option = 0;  % Discretization Option: whether to use c2d (1) or derive manually (0)
ode45_option = 0;           % Simulation Option: whether to use ode45 (1) or manual RK4 (0)

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
Cc =  [0, 0, 1, 0];
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

%% Optimal Feedback Gain Computation
G1 = [Ad,      zeros(nx,nr);
      -Cd*Ad, eye(nr)];
H1 = [Bud; -Cd*Bud];

[P,KK,~] = idare(G1,H1,Q,R,[],[]);
K = KK(1,1:end-nr);
Ki = -KK(1,end-nr+1: end);

%% Main simulation loop
Xlog = nan(nx, Nsim+1);
Xlog(:,1) = x0;
Ucomplog = nan(nu, Nsim);
Uactlog = nan(nu, Nsim);
Vlog = zeros(nr, Nsim+1);
dFlog = zeros(1,Nsim);
reflog = nan(1, Nsim);

for k = 1:Nsim
    % % Reference Update
    if k > (10/Ts)
        xr = -0.5;
    end
    reflog(:,k) = xr;
    % % Disturbance Update
    % % (k >= (20/Ts)) && (k <= (20/Ts)+1) % > (20/Ts) && k <= (25/Ts)
    % if ((k-1)*Ts > 20) &&((k-1)*Ts <= 20.5)
    %     dk = 2*((k-1)*Ts) - 40;
    % elseif ((k-1)*Ts > 20.5) &&((k-1)*Ts <= 21)
    %     dk = 1;
    % elseif ((k-1)*Ts > 21) &&((k-1)*Ts <= 21.5)
    %     dk = -2*((k-1)*Ts) + 43;
    % else
    %     dk = 0;
    % end

    if (k >= (20/Ts)) && (k <= (20/Ts)+10)
        dk = 1;
    else
        dk = 0;
    end

    dFlog(:,k) = dk;

    % % Current Output
    yk = Cd*Xlog(:,k);
    % % Integral State
    Vlog(:,k+1) = Vlog(:,k) + xr - yk;
    
    % Computed Control (unbounded)
    Ucomplog(:,k) = -K*Xlog(:,k) + Ki*Vlog(:,k+1);
    % Actual Control (saturated)
    Uactlog(:,k) = max(u_min,min(Ucomplog(:,k),u_max));
    % System Update
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
legend({'LQR-i','Constraints'},'Interpreter','latex');
ylabel('$\theta (\mathrm{rad})$','Interpreter','latex');
title('(a) Rod Angle','Interpreter','latex')
% ylim([-0.6 0.40]); yticks([-0.55 -0.2618 0 0.2618]);
box on; grid on;
ax = gca; ax.FontSize = 14;

nexttile
hold on; box on; grid on;
plot(Tspan, Xlog(3,1:end-1),'LineStyle','-','LineWidth',lw,'Color','b');
plot(Tspan,reflog,'LineStyle','--','LineWidth',lw,'Color','k');
legend({'LQR-i','Reference'},'Interpreter','latex');
ylabel('$x (\mathrm{m})$','Interpreter','latex');
xlabel('Time (s)','Interpreter','latex');
title('(b) Cart Position','Interpreter','latex')
% ylim([-1 1.2]); yticks([-0.9 -0.5 0 1]);
ax = gca; ax.FontSize = 14;

fig1.TileSpacing = 'compact';
fig1.Padding = 'compact';

% % % % ------------------------------------------------- % % % %

figure(2)
fig2 = tiledlayout(2,1);

nexttile
stairs(Tspan,Vlog(1,1:end-1),'LineStyle','-','LineWidth',lw,'Color','b');
ylabel('$v (\mathrm{m})$','Interpreter','latex');
title('(a) Integral State','Interpreter','latex');
% ylim([-90 120]); yticks([-54 0 107]);
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
% ylim([-34 52]); yticks([-30 0 30]);
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