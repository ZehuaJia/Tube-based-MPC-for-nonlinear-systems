function [P, K, lambda, mu] = NMPC_get_RIS (A, B, Bw)
%% Get the RIS 
% By Zehua Jia, jiazehua@sjtu.edu.cn

%% System description
% x' = A * x + g(x) + B * u + Bw * w
% g(x) = [0; -0.25 * x2^3];
% The designed local state-feeback control law is u = K * x.

%% Some Tips
% A. The constraints should be Linear, which means nonlinear terms are not
% permitted (such as x*y and x^2, where x,y are both sdpvars). 

%% Solve the optimization problem using YALMIP
% Define decision matrix variables
[n,~] = size(A);
[~,m] = size(B);
mu = sdpvar(1);
X = diag(sdpvar(n,1)); 
Y = sdpvar(m,n,'full'); 
% lambda0 = sdpvar(1); 
lambda0 = 1.8;

% Define LMI constraints
LMI = [(A * X + B * Y)' + A * X + B * Y + lambda0 * X, Bw; Bw', -mu];
% const = [LMI <= 0, X >= 0, mu >=0, lambda0 >= 1.8];
const = [LMI <= 0, X >= 0, mu >=0];


%% Solving
% optimize(const,obj,sdpsettings('solver','Mosek'))
solvesdp(const)
%% Get results
X1 = value(X);
Y1 = value(Y);
P = X1^(-1);
K = Y1 * P;
% lambda = value(lambda0);
lambda = lambda0;
mu = value(mu);
end
