% Noisy Linear System (NsyLS)
% Part I - KF + OMP, Part II - KF + POMP
clear; clc; close all
rng(0)
% Simulation Parameters
NoisedB = [-10,-20,-30,-40]; sig_v = 10.^(NoisedB/20);
lns = length(NoisedB);
NSys = 100; n = 80; m = 160; p = n; sig_w = sig_v;
S = floor(5:10:n);
ls = length(S);
K = n/2;

% Initialization
Xf = 10*rand(n,NSys);
x0 = zeros(n,1); R_x = sqrt(1)*eye(n);% Reachability
initPertb = randn(n,NSys);
NMSEi1 = zeros(ls,lns,NSys); % NMSE for OMP

% System Matrices
MA = Erdos_Renyi(n, NSys);
MB = randn(n,m,NSys); MC = rand(p,n,NSys);

tic;
parfor i=1:NSys
    u_omp = zeros(m,1);
    A = MA(:,:,i); B = MB(:,:,i); C = MC(:,:,i);
    xf = Xf(:,i);
    
    for ns=1:lns %index for noise levels
        v_sd = sig_v(ns); w_sd = v_sd;
        V = (v_sd^2)*eye(n); W = (w_sd^2)*eye(p);
        R_v = v_sd*eye(n); R_w = w_sd*eye(p);
        v = R_v*randn(n,K); w = R_w*randn(p,K);
        
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
                X1(:,j+1) = A*X1(:,j) + B*u_omp + v(:,j);
                y1 = C*X1(:,j) + w(:,j);
                
                % New State Prediction
                [X_prd1(:,j+1), P1] = KF_prd(A,B,C,V,W,P1,X_prd1(:,j),y1,u_omp);
                %[X_prd2(:,j+1), P2] = KF_prd(A,B,C,V,W,P2,X_prd2(:,j),y2,u_pomp2);
                %}

                %
                % -- Estimation based Control
                w = randn(p,1);
                y1 = C*X1(:,j) + R_w*w;
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
            NMSEi1(l,ns,i) = vecnorm(xf-X1(:,end)).^2;
        end
    end
end
toc;

NMSE1 = 10*log10(sum(NMSEi1,3)/NSys);
%% Plotting MSE vs Sparsity -Varying Noise Lvls
figure();
plot(S,NMSE1.','LineWidth',3);
grid on;
str = sprintf('\\sigma^2 = ');
legend(strcat(str,num2str(NoisedB.'),' dB'));
ylabel("OMP $10\log{\bf E ||x_f-x||_2^2}$",'Interpreter','latex');
xlabel('Sparsity S');
set(gca,'FontSize',20,'FontWeight','bold');
str = sprintf('OMP n=%d, m=%d, p=%d, NSys = %d',n,m,p,NSys);
title(str,'FontSize',15,'FontWeight','normal');
xlim([0, n]);