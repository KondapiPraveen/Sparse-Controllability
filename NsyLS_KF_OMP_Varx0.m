% Noisy System - MSE vs ||x_0||^2
clear; clc; close all
rng(0)
NTr = 100; n = 80; m = 2*n; p = n; %sig_v = 1e-3; sig_w = sig_v;
NoisedB = [-10,-20,-30,-40]; sig_v = 10.^(NoisedB/20);
lns = length(NoisedB);
s = floor((0.05)*n);
%S = floor(10:5:n); %ls = length(S);
NoisedB = 20*log10(sig_v);
fct = [0,3,6,9,12,15,18,21]; lfc = length(fct);

% Covariance Matrices
%V = (sig_v^2)*eye(n); W = (sig_w^2)*eye(p);
% R_v = sig_v*eye(n); R_w = sig_w*eye(p);
R_x = sqrt(1)*eye(n);
K = n/2; % # Time Steps

NMSEi1 = zeros(lns,lfc,NTr);
X0 = randn(n,NTr); X0 = X0./vecnorm(X0); % Initial State
xf = 1*ones(n,1); % Reachability
initPertb = zeros(n,NTr); % initPertb = randn(n,NTr);

A = Erdos_Renyi(n,1); B = randn(n,m); C = rand(p,n);

tic;
parfor i=1:NTr % Iteration over # Trials
    for f=1:lfc % Norm of x0
        x0 = f*X0(:,i);
        for ns=1:lns % various Noise Levels
            v_sd = sig_v(ns); w_sd = v_sd;
            V = (v_sd^2)*eye(n); W = (w_sd^2)*eye(p);
            R_v = v_sd*eye(n); R_w = w_sd*eye(p);
            v = R_v*randn(n,K); w = R_w*randn(p,K);
        
            X1 = zeros(n,K+1); % Record of the State Vectors (OMP)
            X_est1 = zeros(n,K); % Record of the Estimated State Vectors (OMP)
            X_prd1 = zeros(n,K); % Record of the Predicted State Vectors (OMP)

            X1(:,1) = x0+R_x*(initPertb(:,i));
            X_prd1(:,1) = x0; P1 = R_x;
            
            for j=1:K % Iteratioin over Time Steps
                %{          
                % -- To Shut the Prediction based Control
                % OMP Input Generation
                x_hat1 = xf-A*X_prd1(:,j);
                u_omp = OMP(B,x_hat1,s);

                % System Update
                v = randn(n,1); w = randn(p,1);
                X1(:,j+1) = A*X1(:,j) + B*u_omp + R_v*v;
                y1 = C*X1(:,j) + R_w*w;

                % New State Prediction
                [X_prd1(:,j+1), P1] = KF_prd(A,B,C,V,W,P1,X_prd1(:,j),y1,u_omp);
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
                u_omp = OMP(B,x_hat1,s);
                %Eu_omp = norm(u_omp)^2;
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
            NMSEi1(ns,f,i) = vecnorm(xf-X1(:,end)).^2;
        end
    end
end
toc;
%
NMSE1 = 10*log10(sum(NMSEi1,3)/NTr);
%% Plotting - MSE vs ||xf|| - Varying Sparsity
figure();
plot(fct,NMSE1.','LineWidth',3);
grid on;
str = sprintf('\\sigma^2 = ');
legend(strcat(str,num2str((sig_v.^2).','%.e')),'NumColumns',2);
ylabel("OMP $10\log{\bf E ||x_f-x||_2^2}$",'Interpreter','latex');
xlabel('$\bf ||x_0||_2$','Interpreter','latex');
set(gca,'FontSize',20,'FontWeight','bold');
str = sprintf('OMP n=%d, m=%d, p=%d, NTr = %d, S=%d',n,m,p,NTr,s);
title(str,'FontSize',15,'FontWeight','normal');
xlim([0, fct(end)]); ylim([-25, 15]);