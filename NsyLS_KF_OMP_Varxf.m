% Noisy System - MSE vs ||x_f||^2
clear; clc; close all
rng(0)
NTr = 100; n = 80; m = 2*n; p = n; sig_v = 1e-2; sig_w = sig_v;
S = floor((0.1:0.1:0.5)*n);
%S = floor(10:5:n);
NoisedB = 20*log10(sig_v);
ls = length(S);
fct = [0,3,6,9,12,15,18,21]; lfc = length(fct);
% CoSaMP parameters
max_itr = 20; tol=0.1;

% Covariance Matrices
V = (sig_v^2)*eye(n); W = (sig_w^2)*eye(p);
R_x = sqrt(1)*eye(n); R_v = sig_v*eye(n); R_w = sig_w*eye(p);
K = n/2; % # Time Steps

NMSEi1 = zeros(ls,lfc,NTr);
Xf = randn(n,NTr); Xf = Xf./vecnorm(Xf); % Final State
x0 = zeros(n,1); % Reachability
initPertb = zeros(n,NTr); % initPertb = randn(n,NTr);

MA = Erdos_Renyi(n,NTr); MB = randn(n,m,NTr); MC = rand(p,n,NTr);

tic;
parfor i=1:NTr
    u_omp = zeros(m,1);
    A = MA(:,:,i); B = MB(:,:,i); C = MC(:,:,i);
    for f=1:lfc
        xf = f*Xf(:,i);
        for l=1:ls
            s = S(l);
            X1 = zeros(n,K+1); % Record of the State Vectors (OMP)
            X_est1 = zeros(n,K); % Record of the Estimated State Vectors (OMP)
            X_prd1 = zeros(n,K); % Record of the Predicted State Vectors (OMP)

            X1(:,1) = x0+R_x*(initPertb(:,i));
            X_prd1(:,1) = x0; P1 = R_x;
            
            for j=1:K
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
                % u_omp = OMP(B,x_hat1,s);
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
                %}
            end
            NMSEi1(l,f,i) = vecnorm(xf-X1(:,end)).^2;
        end
    end
end
toc;
%
NMSE1 = 10*log10(sum(NMSEi1,3)/NTr);
%% Plotting - MSE vs ||xf|| - Varying Sparsity
figure();
plot(fct,NMSE1.','LineWidth',3,'MarkerSize',10);
grid on;
legend(strcat('$s$ = ',num2str(S.')),'NumColumns',2,'Interpreter','latex');
ylabel("OMP $10\log{\bf E ||x_f-x||_2^2}$",'Interpreter','latex');
xlabel('$\bf ||x_f||_2$','Interpreter','latex');
set(gca,'FontSize',20,'FontWeight','bold');
str = sprintf('OMP n=%d, m=%d, p=%d, NTr = %d, \\sigma^2=%d dB',n,m,p,NTr,NoisedB);
title(str,'FontSize',15,'FontWeight','normal');
xlim([0, fct(end)]);