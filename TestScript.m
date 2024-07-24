%% POMP Script Test
% Not a great Recovery but reduces the residual
M = 8; N = 256; L = 8; % L sparse vectors are concatanated.
D = randn(M,N); D = D./vecnorm(D); s = zeros(N,1);
Ki = 5; K = L*Ki;
I = (0:L-1)*(N/L);
Supp = unique(ceil((N/L)*(rand(Ki,L)+[0:L-1])));
s(Supp(:)) = randn(numel(Supp),1);
b = D*s;
sk = POMP(D,b,L,Ki);
norm(s-sk)

%% SpaIpDsg (Sparse Input Design) Script Test
clear ;clc; close all;
n = 5; m = 5;
A = randn(n); B = randn(n,m);
q = length(minpoly(A))-1;
S_min = n-rank(A)+1;
s = max(S_min, ceil(m/5));
R_Bs = min(rank(B),s);
%
x0 = rand(n,1); xf = rand(n,1);
R = CtrlMatrix(A,B,n);
[~, ~, Residue] = SpaIpDsg(x0,xf,A,B,s,q);
Kst = ceil(n/R_Bs);
semilogy(Kst+(0:length(Residue)-1), Residue, 'LineWidth', 1.5, 'Marker', 's', 'MarkerSize', 10)
title(['Sparse controllability with s = ', num2str(s)])
xlabel('Time Steps')
ylabel('Residue State Error')

%% Greedy Scheduling Test
clear ;clc; close all;
n = 10; m = 10;% # States and # Inputs
A = randn(n); B = randn(n,m);
s = ceil(m/2);
ts = n; % # Time Steps
[S_G, itr] = GreedyScheduling_Aopt(A,B,ts,s);
% Random Sampling Test (Run along with Greedy Scheduling)
[S_R, W] = RandSamp_Aopt(A,B,ts,s);
R = CtrlMatrix(A,B,ts);
W_S = R(:,S_R)*(R(:,S_R).');
% fprintf('The rank of Sparse Scheduled Controllability Gramian is %d (Random Sampling) \n',rank(W_S))
% fprintf('The Condition no. of Controllability Matrix: %d (Random Sampling) \n',cond(W_S))
%% Kalman Filter Test
n = 10; m = n; p = 15; TSteps = 10; sig_v = 1e-2; sig_w = 1e-2;
V = (sig_v^2)*eye(n); W = (sig_w^2)*eye(p);
A = randn(n); B = rand(n); C = randn(p,n);
x = zeros(n,TSteps+1); x_est = zeros(n,TSteps); u = rand(m,TSteps);
x(:,1) = randn(n,1); x_0 = zeros(n,1);
for i=1:TSteps
    y = C*x(:,i) + sig_w*randn(p,1);
    x(:,i+1) = A*x(:,i) + B*u(:,i) + sig_v*randn(n,1);
    if i==1
        P_prd = eye(n);
        K = P_prd*(C.')/(C*P_prd*(C.')+W);
        x_est(:,1) = x_0 + K*(y-C*x_0); % State Estimate
        P = (eye(n)-K*C)*P_prd; % Posterior Covariance
    else
        [x_est(:,i),P] = KF(A,B,C,V,W,P,x_est(:,i-1),y,u(:,i-1));
    end
end
figure();
plot(1:TSteps, vecnorm(x(:,1:end-1)-x_est),'red','LineWidth',3)
figure();
plot(1:TSteps,x(1,1:end-1),'red','LineWidth',3,'DisplayName','True')
hold on; grid on
plot(1:TSteps,x_est(1,:),'blue--','LineWidth',3,'DisplayName','Estimate')
legend();
%% OFFSet Evaluation in Noiseless System
clear ; clc; close all;
rng(0);
n = 40; m = 40; K = n; NSys = 1; btch = 15; % State, Input Dimension, Initial Control Time
x0 = zeros(n,1); xf = 10*rand(n,1); % Initial and Final States
MA = Erdos_Renyi(n,NSys); MB = rand(n,m,NSys);
N = 10; % Control Horizon after reaching xf
RegK = zeros(K*n,K*n);
RegN = zeros(N*n,N*n);
S = (0.1:0.1:0.3)*m;
ls = length(S);
u = zeros((K+N*btch)*m,ls,NSys);
x = zeros((K+N*btch)*m,ls,NSys);
Eu = zeros(N*btch+K,ls,NSys); % Control Energy Profile
EX = zeros(N*btch+K,ls,NSys); % State Error Profile
tic;
for k=1:NSys
    A = MA(:,:,k); B = MB(:,:,k);
    AK = A^K; b = xf - AK*x0;
    [M1K, M2K] = bck_lwr_mtx(A,B,RegK,K);
    [M1N, M2N] = bck_lwr_mtx(A,B,RegN,N);
    M2K = M2K(1:K*n,:); M2N = M2N(1:N*n,:);
    R = CtrlMatrix(A,B,K);
    for i=1:ls
        s = S(i);
        u(1:K*m,i,k) = POMP(R,b,K,s);
        x(1:K*n,i,k) = M1K*x0 + M2K*u(1:K*m,i,k);
        x02 = x((K-1)*n+1:K*n,i,k);
        b2 = repmat(xf,N,1) - M1N*x02;
        for j=1:btch
            u((K+N*(j-1))*m+1:(K+N*j)*m,i,k) = POMP(M2N,b2,N,s);
            x((K+N*(j-1))*n+1:(K+N*j)*n,i,k) = M1N*x02 + M2N*u((K+N*(j-1))*m+1:(K+N*j)*m,i,k);
            x02 = x((K+N*j-1)*n+1:(K+N*j)*n,i,k);
            b2 = repmat(xf,N,1) - M1N*x02;
        end
        U = reshape(u(:,i,k),m,N*btch+K);
        Eu(:,i,k) = vecnorm(U).';
        X = reshape(x(:,i,k),n,N*btch+K);
        EX(:,i,k) = vecnorm(xf-X).';
    end
end
toc;

% Plotting
AEX = 10*log10(mean(EX.^2,3)); AEu = mean(Eu.^2,3);
figure()
t  = 1:(N*btch+K);
plot(t,AEu,'LineWidth',2)
legend(strcat('s = ',num2str(S.')));
xlabel('time (T)'); ylabel('$10\log{\bf  E||u||^2}$','Interpreter','latex');
str = sprintf('n = %d, m = %d, K = %d, N = %d, NSys = %d',n,m,K,N,NSys);
title(str);

figure()
plot(t,AEX,'LineWidth',2)
legend(strcat('s = ',num2str(S.')));
xlabel('time (T)'); ylabel('$10\log{\bf E||xf-x||^2}$','Interpreter','latex');
title(str); ylim([0 40])
%% Test FullBI.m (Find LI Columns that minimize Tr(W_S^{-1}))
clc;
n=100; rng(1);
A = zeros(n); r = 15;
for i =1:100
    A2 = Erdos_Renyi(n-15,1); A(1:n-15,1:n-15) = A2;  B = rand(n);
    s = max(r,n-rank(A)); K = ceil(n/s);
    [S,S_k]=FullBLI(A,B,K,s);
end