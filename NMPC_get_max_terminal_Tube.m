function [P, K, alpha] = NMPC_get_max_terminal_Tube (Q, R, ucon, xcon, alpha)
%% Get the maximal terminal region 
% Method adopted from Wenhua Chen et al. 2003 in IJACSP. By Zehua Jia, jiazehua@sjtu.edu.cn

%% Input
% ucon: |u| <= ucon; 
% xcon: |x| <= xcon
% Here x(i) has the same upper and lower bounds, i = 1, 2.
% alpha: x' * P * x <= alpha
%% System description
% x' = A * x + g(x) + B * u + Bw * w
% g(x) = [0; -0.25 * x2^3];
% The designed local state-feeback control law is u = K * x.

%% Some Tips
% A. The constraints should be Linear, which means nonlinear terms are not
% permitted (such as x*y and x^2, where x,y are both sdpvars).

% B. The inverse of a sdpvar is not desired in constraints or objective.
% logdet(X^(-1)) = -logdet(X).

% C. logdet is concave, and then -logdet is convex.

% D. One should always try to reformulate the problem to obtain a convex 
% problem. 

% E. -logdet(a * X) (sdpvar a, X = sdpvar(n,n)) can be reformulated as 
% -logdet(Y) (sdpvar a, Y = sdpvar(n,n)). Then multiply LMI constraints by 
% a in both sides, and replace Y = a * X in constraints, 

%% Initiallization
umax = ucon;
umin = -ucon;
xmax = xcon;
xmin = -xcon;
%% LDI approximation within the selected set x <= 2
N = 2; % The number of LDI 
A = [-1, 2; -3, 4];
B = [0.5; -2];
g(:, :, 1) = [0, 0; 0, -0.75 * xmax^2];
g(:, :, 2) = [0, 0; 0, 0];
for i = 1:N
    F(:,:,i) = [A + g(:, :, i), B];
end

%% Solve the optimization problem using YALMIP
% Express the state/Input constraints in standard form
c(1, :) = [1/xmax, 0];% for x1 < 2
c(2, :) = [0, 1/xmax];% for x2 < 2
c(3, :) = -c(1, :);% for x1 > -2
c(4, :) = -c(2, :);% for x2 > -2
c(5, :) = zeros(1, 2);
c(6, :) = zeros(1, 2);
d(1) = 0;
d(2) = 0;
d(3) = 0;
d(4) = 0;
d(5) = 1/umax;
d(6) = - d(5);
[~, Nc] = size(d); % Number of constraints

% Define decision matrix variables
[n,~] = size(Q);
[~,m] = size(R);
alpha0 = sdpvar(1);
W1 = sdpvar(n,n); % W1 = alpha * W1 (The latter W1 is the W1 in the paper)
W2 = sdpvar(m,n,'full'); % W2 = alpha * W2 (The latter W2 is the W2 in the paper)
W = [W1, W2'];
MAT1 = sdpvar(2*n+m, 2*n+m, N);
MAT2 = sdpvar(1+n, 1+n, Nc);

% Define LMI constraints
for i = 1:N
    MAT1(:,:,i) = [-F(:,:,i) * W' - W * F(:,:,i)', [W1 * Q^(0.5), W2'];
        [W1 * Q^(0.5), W2']',alpha0 * [eye(2), zeros(2,1); zeros(1,2), R^(-1)]]; % Multiplying alpha in both sides in LMI (19)
end

for i = 1:Nc
    MAT2(:,:,i) = [1, c(i, :) * W1 + d(i) * W2; (c(i, :) * W1 + d(i) * W2)', W1]; % Multiplying alpha in both sides in LMI (20)
end

const = [];
for i = 1:N
    const = [const, MAT1(:,:,i)>=0];
end

for i = 1:Nc
    const = [const, MAT2(:,:,i)>=0];
end
% const = [const, alpha <= 10]; % Without this constraint, the solver will turn out error "lack of progress"
obj = -logdet(W1);
const1 = [const, alpha0 <= alpha]; % The LMI method with limited alpha
%% Solving
if alpha == -1
    optimize(const, obj, sdpsettings('solver','sdpt3')) % When input alpha == -1, it means no selected alpha
elseif alpha <= 0
    error('The alphaM must be larger than 0')
else
    optimize(const1, obj, sdpsettings('solver','sdpt3'))
end
%% Get results
W1 = double(W1) / double(alpha0);
W2 = double(W2) / double(alpha0);
P = W1^(-1);
K = W2 * W1^(-1);
alpha = double(alpha0);

%% Plot the obtained ellipsoid (terminal region)
% draw_ellip(P, alpha, 'k')
% hold on % Used for multiple eliipses comparison
end
