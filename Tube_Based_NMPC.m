%% Tube MPC scheme based on robust control invariant set with application to Lipschitz nonlinear systems (2013)
% Tube-based NMPC by Zehua Jia, jiazehua@sjtu.edu.cn

% Simulation reproduction of Yu et. al 2013.

%% System description
% x' = A * x + g(x) + B * u + Bw * w

%% Initiallization
clear all
close all
clc

% rng(0); % fix random seed
%% Model description
A= [-1, 2; -3, 4];
B = [0.5; -2];
Bw = [0; 1];
% g(x) = [0; 0.25 * x2^3];
%% Disturbance
wmax = 0.1;
W_vertex = [-wmax, 0; wmax, 0];
W = Polyhedron(W_vertex);

%% Problem formulation
dt = 0.1; % The sampling period
N = 15; % The prediction horizon N * dt
T = 100; % The simulation time T * dt

Q = diag([0.5, 0.5]);
R = 1;

%% Initial settings
x0 = [3.5; -2.5]; % x0
uM = 2; % The constraint of u for set M of the LDI
xM = 2; % The constraint of x for set M of the LDI
alphaM = 10; % alphaM is the given upper bound of alpha. -1 means no given bound.
% Define the variables storing the actual states and controls
xc = zeros(2, T);
uc = zeros(1, T);
xn = zeros(2, T); % store nominal states
un = zeros(1, T);
tilde_x = zeros(2, T);
xc(:, 1) = x0;
xn(:, 1) = x0;
w = zeros(1, T);

%% Caculate the robust invariant set (RIS)
[P, K, lambda0, mu] = NMPC_get_RIS (A, B, Bw); 
% lambda = lambda0 - 1; 
K1 = [-1.3693, 5.1273]; % The result in the paper
P1 = diag([39.0251, 486.0402]); % The result in the paper
draw_ellip(P, mu * wmax^2/lambda0, 'g') % draw the ellipse
hold on
draw_ellip(P1, 1, 'r')
% Omega = Polyhedron([-0.1127,-0.04581; -0.1127,0.04581; 0.1126,-0.04581; 0.1126,0.04581]);
%% Calculate the Minkowski difference set
x = sdpvar(2,1);
const = x' * P * x <= (mu * wmax^2 / lambda0);
% const = x' * P1 * x <= 1;
obj1 = K * x;
% obj2 = [1, 0] * x;
% obj3 = [0, 1] * x;
optimize(const, obj1);
u_0 = value(obj1);
% optimize(const, obj2);
% x10 = value(obj2);
% optimize(const, obj3);
% x20 = value(obj3);
%% Get the terminal terminal region Xf and terminal penalty for nominal MPC 
[P_nom, K0_nom, alpha_nom] = NMPC_get_max_terminal_Tube (Q, R, uM + u_0, xM, alphaM); 
P_nom1 = [7.9997, -12.2019; -12.2019, 27.0777];
% pause(1)
draw_ellip(P_nom, alpha_nom, 'g.')
hold on
draw_ellip(P_nom1, alpha_nom, 'r.')
%% MPC Optimization problem using YALMIP
for v = 1 : 6 % The loop for collecting data under different disturbances
    x0 = [3.5; -2.5]; % x0
    xc = zeros(2, T);
    uc = zeros(1, T);
    xn = zeros(2, T); % store nominal states
    un = zeros(1, T);
    tilde_x = zeros(2, T);
    xc(:, 1) = x0;
    xn(:, 1) = x0;
    w = zeros(1, T);
    for i = 1 : T-(N-1)
        % Define decision variables
        x = sdpvar(2, N);
        u = sdpvar(1, N-1);
        w(i) = pick_random_disturbance(W);
%             w(i) = 0.1;
        % Define constraints
        const = [u <= (uM + u_0), u >= (-uM - u_0), x(:,1) == x0];
        %     const = [u <= 2, u >= -2, x(:,1) == x0];
        for k = 1 : N-1
            const = [const, x(:,k+1) == x(:, k) + dt * ( A * x(:, k) + B * u(k) + [0; -0.25 * x(2, k)^3] )];
        end
        const = [const, x(:, N)' * P_nom * x(:, N) <= alpha_nom];
        % Define objective
        obj = 0;
        for j = 1 : N-1
            obj = obj + x(:, j)' * Q * x(:, j) + u(j)' * R * u(j);
        end
        obj = obj + x(:, N)' * P_nom * x(:, N);
        % Optimization
        optimize(const, obj);
        %     optimize(const, obj, sdpsettings('solver','fmincon','fmincon.maxiter',3000));
        % Control and updates
        un(i) = value(u(1));
        uc(i) = un(i) + K * tilde_x(:, i);
        xn(:, i+1) = xn(:, i) + dt * ( A * xn(:, i) + B * un(i) + [0; -0.25 * xn(2, i)^3] );
        xc(:, i+1) = xc(:, i) + dt * ( A * xc(:, i) + B * uc(i) + [0; -0.25 * xc(2, i)^3] + Bw * w(i));
        tilde_x(:, i+1) = xc(:, i+1) - xn(:, i+1);
        x0 = xn(:, i+1);
        %      x0 = value(x(:,2));
        i % Show the current step
        if v == 1
            figure(4)
            draw_ellip2(P, mu * wmax^2/lambda0, [xn(1,i), xn(2,i)] ,'y')
            hold on
        end
    end
    
    for i = T-N+2 : T-1
        w(i) = pick_random_disturbance(W);
        % w(i) = 0.1;
        tilde_x(:, i) = xc(:, i) - xn(:, i);
        un(i) = value(u(i-T+N));
        xn(:, i+1) = value(x(:,i-T+N));
        uc(i) = un(i) + K * tilde_x(:, i);
        xc(:, i+1) = xc(:, i) + dt * ( A * xc(:, i) + B * uc(i) + [0; -0.25 * xc(2, i)^3] + Bw * w(i));
    end
    tilde_x(:, T) = xc(:, T) - xn(:, T);
    
    % Save the data under different disturbances
    % pathname = 'F:\สตั้สา\ทยีๆ\MPC\Reproduction of TMPC\tube-based test data set\';
    % str0 = num2str(v);
    % filename1 = ['xc',str0];
    % filename2 = ['xn',str0];
    % filename3 = ['uc',str0];
    % filename4 = ['un',str0];
    % save([pathname,filename1],'xc');
    % save([pathname,filename2],'xn');
    % save([pathname,filename3],'uc');
    % save([pathname,filename4],'un');
    t = 0: dt: (T-1) * dt;
    t_com = 0 : dt: (N-1) * dt;
    figure(2)
    plot(t, xc(1,:), 'k', t, xc(2,:), 'b')
    hold on
    plot(t, xn(1,:), 'r', t, xn(2,:), 'g')
    hold on
    
    
    figure(3)
    plot(t, uc(:))
    hold on
    plot(t, un(:),'r')
    hold on
    
    figure(4)
    draw_ellip(P_nom, alpha_nom, 'g')
    hold on
    % plot(Omega)
    % hold on
    % draw_ellip(P, mu * 0.01/lambda, 'y')
    % hold on
    plot(xn(1, :),xn(2, :),'-*r')
    hold on
    plot(xc(1, :),xc(2, :),'-*')
    hold on
    % xlim([-4,4])
    % ylim([-5,2])
    
    figure(5)
    plot(t, tilde_x(1,:), 'k', t, tilde_x(2,:), 'b')
    hold on
    
    
    figure(6)
    plot(tilde_x(1,:), tilde_x(2,:), 'b*')
    hold on
    draw_ellip2(P, mu * wmax^2/lambda0, [0, 0] ,'r.')
end

%% Plot the trajectory
% t = 0: dt: (T-1) * dt;
% t_com = 0 : dt: (N-1) * dt;
% figure(1)
% plot(t, xc(1,:), 'k', t, xc(2,:), 'b')
% hold on
% plot(t, xn(1,:), 'r', t, xn(2,:), 'g')
% 
% 
% figure(2)
% plot(t, uc(:))
% hold on
% plot(t, un(:),'r')
% 
% figure(3)
% draw_ellip(P_nom, alpha_nom, 'g')
% hold on
% % plot(Omega)
% % hold on
% % draw_ellip(P, mu * 0.01/lambda, 'y')
% % hold on
% plot(xn(1, :),xn(2, :),'-*r')
% hold on
% plot(xc(1, :),xc(2, :),'-*')
% % xlim([-4,4])
% % ylim([-5,2])
% 
% figure(4)
% plot(t, tilde_x(1,:), 'k', t, tilde_x(2,:), 'b')
% hold on
% 
% 
% figure(5)
% plot(tilde_x(1,:), tilde_x(2,:), 'b*')
% hold on
% draw_ellip2(P, mu * wmax^2/lambda, [0, 0] ,'r.')




