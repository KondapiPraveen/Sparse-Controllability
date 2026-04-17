% Fixed Point Algorithms for the Ricatti Recursion

% Parameters
A = eye(3); % Transfer matrix
C = eye(3); % Output matrix
V = sqrt(1)*eye(3); % Process noise covariance
W = sqrt(1)*eye(3);% Measurement noise covariance
n = size(A,1); tol = 1e-3;

% Compute k_1,k_2,k_3 with s=1
% k1 computation
eta = 1e-1; d0 = 30;
zi = 1-1/n; 
k1 = ceil(log(eta/(3*d0^2))/log(zi));
eps = eta*(1-zi)/(3*(1+3*zi));

% Init
R_post_new = zeros(n,n);

canItr = 1; count=0;

while canItr
    R_pre = A*R_post_new*A.' + V;
    L = R_pre*C.'/(C*R_pre*C.' + W);
    R_post_old = R_post_new;
    R_post_new  = (eye(n)-L*C)*R_pre;
    canItr = norm(R_post_new - R_post_old,'fro') > tol;
    count = count+1;
end

% k2 computation
P = R_post_new;

canItr = 1; k2 = -1;
R_pre = zeros(n,n);
R_post_old = zeros(n,n);
R_post_new = zeros(n,n);
cst_k3 = 0;
while canItr
    R_pre = A*R_post_new*A.' + V;
    L = R_pre*C.'/(C*R_pre*C.' + W);
    R_post_old = R_post_new;
    R_post_new  = (eye(n)-L*C)*R_pre;
    canItr = abs(trace(R_post_new) - trace(P)) > eps;
    cst_k3 = cst_k3 + trace(R_post_new);
    k2 = k2+1;
end

k3 = k2 + ceil(log(eps/cst_k3)/log(zi));