% Noisy Linear System (NsyLS) - MSE vs Time Steps(k)
% Part I - KF + OMP, Part II - KF + POMP
clear; clc; close all
rng(0)
% Simulation Parameters
NSys = 100; n = 80; m = 2*n; p = n; sig_v = 1e-2; sig_w = sig_v;
%S = floor(([0.1:0.05:0.25, 1])*m);
S = ceil((0.1:0.1:0.5)*n);
NoisedB = 20*log10(sig_v);
%S = m;
ls = length(S);
% CoSaMP parameters
max_itr = 20; tol=0.1;

% Covariance Matrices
V = (sig_v^2)*eye(n); W = (sig_w^2)*eye(p);
R_x = sqrt(1)*eye(n); R_v = sig_v*eye(n); R_w = sig_w*eye(p);
K = ceil(n/2); % # Time Steps
k = 5; % Control Horizon

% Initialization
Xf = 10*rand(n,NSys); % x0 = Xf;
x0 = zeros(n,1); initPertb = zeros(n,NSys); %initPertb = randn(n, NSys); % Reachability
%Xf = zeros(n,1); x0 = 10*rand(n,NSys); % Controllability
NMSEi1 = zeros(ls,K+1,NSys); NMSEi2 = zeros(ls,K+1,NSys); % OMP, POMP
UOMP = zeros(ls,K,NSys); UPOMP = zeros(ls,K,NSys);
%R = zeros(m*k); % xf = zeros(2*n,1);
R = 0*eye(m*k); % Regularizer for Input

% System Matrices
MA = Erdos_Renyi(n, NSys);
MB = randn(n,m,NSys); MC = rand(p,n,NSys);

tic;
parfor i=1:NSys
    u_omp = zeros(m,1);
    u_pomp2 = zeros(m,1);
    x_hat2 = zeros((n+m)*k,1); % Initialization;
    A = MA(:,:,i); B = MB(:,:,i); C = MC(:,:,i);
    xf = Xf(:,i); % Desired final State
    rxf = repmat(xf,k,1); % rxf(1:n*k) = rxf1;
    
    % xf(1:n) = xf1;
    
    [M1,M2] = bck_lwr_mtx(A,B,R,k);
    for l=1:ls
        s = S(l);
        X1 = zeros(n,K+1); % Record of the State Vectors (OMP)
        X2 = zeros(n,K+1); % Record of the State Vectors (POMP)
        X_est1 = zeros(n,K); % Record of the Estimated State Vectors (OMP)
        X_est2 = zeros(n,K); % Record of the Estimated State Vectors (POMP)
        X_prd1 = zeros(n,K); % Record of the Predicted State Vectors (OMP)
        X_prd2 = zeros(n,K); % Record of the Predicted State Vectors (POMP)
        
        %Ak = A^k; R = CtrlMatrix(A,B,k);
        
        X1(:,1) = x0+R_x*(initPertb(:,i)); X2(:,1) = X1(:,1); % Reachability
        %X1(:,1) = x0(:,i); X2(:,1) = x0(:,i); % Controllability
        %X_est1(:,1) = x0(:,i); X_est2(:,1) = x0(:,i);
        X_prd1(:,1) = x0; X_prd2(:,1) = x0;
        P1 = R_x; P2 = R_x;
        
        %{
        y = C*X(:,1) + R_w*randn(p,1);
        P_prd = R_x;
        K = P_prd*(C.')/(C*P_prd*(C.')+W); % Kalman Gain
        X_est1(:,1) = x0 + K*(y-C*x0);
        X_est2(:,1) = X_est1(:,1);
        % Calculate Inputs
        x_hat1 = xf-A*X_est1(:,1);
        u_omp = OMP(B,x_hat1,S);
        x_hat2 = xf-Ak*X_est(:,1);
        u_pomp = POMP(R,x_hat2,n,S);
        u_pomp2 = u_pomp(1:m);
        %}
        
        for j=1:K
            %{          
            % -- To Shut the Prediction based Control
            % OMP Input Generation
            x_hat1 = xf-A*X_prd1(:,j);
            u_omp = OMP(B,x_hat1,s);
            Eu_omp = norm(u_omp)^2;
            UOMP(l,j,i) = Eu_omp;
            %{
            % -- To Shut POMP
            % POMP Input Generation
            x_f2 = rxf-M1*X_prd2(:,j);
            x_hat2(1:n*k) = x_f2;
            u_pomp = POMP(M2,x_hat2,k,s);
            u_pomp2 = u_pomp(1:m);
            Eu_pomp2 = norm(u_pomp2)^2;
            UPOMP(l,j,i) = Eu_pomp2;
            %}
            % System Update
            v = randn(n,1); w = randn(p,1);
            X1(:,j+1) = A*X1(:,j) + B*u_omp + R_v*v;
            y1 = C*X1(:,j) + R_w*w;
            %{        
            % -- To Shut POMP
            X2(:,j+1) = A*X2(:,j) + B*u_pomp2 + R_v*v;
            y2 = C*X2(:,j) + R_w*w;
            %}
            % New State Prediction
            [X_prd1(:,j+1), P1] = KF_prd(A,B,C,V,W,P1,X_prd1(:,j),y1,u_omp);
            %[X_prd2(:,j+1), P2] = KF_prd(A,B,C,V,W,P2,X_prd2(:,j),y2,u_pomp2);
            %}
            
            %
            % -- Estimation based Control
            w = randn(p,1);
            y1 = C*X1(:,j) + R_w*w; % OMP
            %y2 = C*X2(:,j) + R_w*w; % POMP
            
            % -- For Controllability shut down the if condition and active
            % else condition for j~=1
            if j==1
                KG = P1*(C.')/(C*P1*(C.')+W); % Kalman Gain
                X_est1(:,1) = X_prd1(:,1) + KG*(y1-C*X_prd1(:,1)); % State Estimate
                P1 = (eye(n)-KG*C)*P1; % Posterior Covariance
                P2 = P1;
            %end
            else
                % New State Estimation
                [X_est1(:,j), P1] = KF(A,B,C,V,W,P1,X_est1(:,j-1),y1,u_omp);
                %[X_est2(:,j), P2] = KF(A,B,C,V,W,P2,X_est2(:,j-1),y2,u_pomp2); % POMP
            end
            
            % OMP Input Generation
            x_hat1 = xf-A*X_est1(:,j);
            %u_omp = OMP(B,x_hat1,s);
            u_omp = CoSaMP(x_hat1,B,s,max_itr,tol);
            Eu_omp = norm(u_omp)^2;
            UOMP(l,j,i) = Eu_omp;
            %{
            % -- To Shut POMP
            % POMP Input Generation
            x_f2 = rxf-M1*X_est2(:,j);
            x_hat2(1:n*k) = x_f2;
            u_pomp = POMP(M2,x_hat2,k,s);
            u_pomp2 = u_pomp(1:m);
            Eu_pomp2 = norm(u_pomp2)^2;
            UPOMP(l,j,i) = Eu_pomp2;
            %}
            % System Update
            v = randn(n,1);
            X1(:,j+1) = A*X1(:,j) + B*u_omp + R_v*v;
            %{          
            % -- To Shut POMP
            X2(:,j+1) = A*X2(:,j) + B*u_pomp2 + R_v*v;
            %}
        end
        NMSEi1(l,:,i) = vecnorm(xf-X1).^2;
        NMSEi2(l,:,i) = vecnorm(xf-X2).^2;
    end
end
toc;

NMSE1 = 10*log10(sum(NMSEi1,3)/NSys); NMSE2 = 10*log10(sum(NMSEi2,3)/NSys);
%% Plotting - OMP vs POMP
%{
figure();
plot(1:K+1,NMSE1,'r','DisplayName','OMP','LineWidth',3)
hold on; grid on;
plot(1:K+1,NMSE2,'b--','DisplayName','POMP','LineWidth',3)
legend();
ylabel("$||x-xf||_2^2$",'Interpreter','latex','FontSize',16);
xlabel('Time Steps K','FontSize',12)
%}
%% Plotting - MSE vs Time Steps - Varying Sparsity
figure();
MkrInd = [1,2,3,4,5,6,7,8,9,11,13,15,17,19,21,23,25];
plot(0:K,NMSE1.','LineWidth',3,'MarkerIndices',MkrInd,'MarkerSize',10);
grid on;
legend(strcat('$s$ = ',num2str(S.')),'Interpreter','latex','NumColumns',2);
ylabel("OMP $10\log{\bf E ||x_f-x||_2^2}$",'Interpreter','latex');
xlabel("Time steps $(k)$",'Interpreter','latex');
set(gca(),'FontSize',20,'FontWeight','bold');
str = sprintf('OMP n=%d, m=%d, p=%d, NSys = %d, \\sigma^2=%d dB',n,m,p,NSys,NoisedB);
title(str,'FontSize',15,'FontWeight','normal');
xlim([0, 25]); xticks([0,5,10,15,20,25])
%% POMP for S
%{
figure();
plot(0:K,NMSE2.','LineWidth',3)
grid on;
legend(strcat('s = ',num2str(S.')));
ylabel("POMP $10\log{\bf E ||xf-x||_2^2}$",'Interpreter','latex');
xlabel('Time Steps K','FontSize',12);
set(gca,'FontSize',20,'FontWeight','bold');
str = sprintf('POMP n=%d, m=%d, p=%d, NSys=%d, N=%d, \\sigma^2=%d dB',n,m,p,NSys,k,NoisedB);
title(str,'FontSize',15,'FontWeight','normal');
xlim([0, K]);
%}
%% ||u||_2^2 for S - OMP
figure();
plot(0:K-1,mean(UOMP,3).','LineWidth',3);
grid on;
legend(strcat('s = ',num2str(S.')));
xlabel('Time Steps K','FontSize',12);
ylabel("OMP $\bf E ||u||_2^2$",'Interpreter','latex');
str = sprintf('OMP n=%d, m=%d, p=%d, NSys = %d, \\sigma^2=%d dB',n,m,p,NSys,NoisedB);
title(str,'FontSize',15,'FontWeight','normal');
xlim([0,K]);
%% ||u||_2^2 for S - POMP
%{
figure();
plot(0:K-1,mean(UPOMP,3).','LineWidth',3);
grid on;
legend(strcat('s = ',num2str(S.')));
xlabel('Time Steps K','FontSize',12);
ylabel("POMP $\bf E ||u||_2^2$",'Interpreter','latex');
str = sprintf('POMP n=%d, m=%d, p=%d, NSys = %d, N=%d, \\sigma^2=%d dB',n,m,p,NSys,k,NoisedB);
title(str,'FontSize',15,'FontWeight','normal');
xlim([0,K]);
%}