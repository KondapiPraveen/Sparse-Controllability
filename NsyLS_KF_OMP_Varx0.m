% Noisy System - MSE vs ||x_0||^2
clear; clc; close all
rng(0)
NTr = 100; n = 80; m = 2*n; p = n; %sig_v = 1e-3; sig_w = sig_v;
NoisedB = [-10,-20,-30,-40]; sig_v = 10.^(NoisedB/20);
lns = length(NoisedB);
S = [4, 10]; ls = length(S);
%s = floor((0.05)*n);
%S = floor(10:5:n);
NoisedB = 20*log10(sig_v);
fct = [0,3,6,9,12,15,18,21]; lfc = length(fct);
% CoSaMP parameters
max_itr = 20; tol=0.1;

% Covariance Matrices
%V = (sig_v^2)*eye(n); W = (sig_w^2)*eye(p);
% R_v = sig_v*eye(n); R_w = sig_w*eye(p);
R_x = sqrt(1)*eye(n);
K = n/2; % # Time Steps

NMSEi1 = zeros(lns,lfc,ls,NTr);
X0 = randn(n,NTr); X0 = X0./vecnorm(X0); % Initial State
xf = 1*ones(n,1); % Reachability
initPertb = zeros(n,NTr); % initPertb = randn(n,NTr);

MA = Erdos_Renyi(n,NTr); MB = randn(n,m,NTr); MC = rand(p,n,NTr);
LwrBnd = zeros(lns, NTr);

tic;
parfor i=1:NTr % Iteration over # Trials
    u_omp = zeros(m,1);
    A = MA(:,:,i); B = MB(:,:,i); C = MC(:,:,i);
    for k = 1:ls
        s = S(k);
        for f=1:lfc % Norm of x0
            x0 = fct(f)*X0(:,i);
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
                    % w = randn(p,1);
                    y1 = C*X1(:,j) + w(:,j); % OMP
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
                    %u_omp = CoSaMP(x_hat1,B,s,max_itr,tol);
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
                    % v = randn(n,1);
                    X1(:,j+1) = A*X1(:,j) + B*u_omp + v(:,j);
                    %{          
                    % -- To Shut POMP
                    X2(:,j+1) = A*X2(:,j) + B*u_pomp2 + R_v*v;
                    %}
                end
                NMSEi1(ns,f,k,i) = vecnorm(xf-X1(:,end)).^2/(norm(xf)^2);
                if f == 1
                    LwrBnd(ns,i) = (trace(V) + trace(A*P1*A.'))/(norm(xf)^2);
                end
            end
        end
    end
    
end
toc;
%
NMSE1 = 10*log10(sum(NMSEi1,3)/NTr);
%% Plotting - MSE vs ||xf|| - Varying Sparsity
DefaultColors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; ...
    0.4660 0.6740 0.1880; 0.3010 0.7350 0.9330; 0.6350 0.0780 0.1840];
figure();
Stylz = ["-","-.",":"];
a1 = axes();
Colors = repmat(DefaultColors(1:lns,:),3,1);
colororder(Colors)
for k=1:ls
    plot(fct, NMSE1(:,:,k).',Stylz(k),'LineWidth',2.5,'MarkerSize',10); hold on
end
plot(fct, ones(numel(fct),1)*lwbnd.','--','LineWidth',2.5,'MarkerSize',10)
grid on;
str = sprintf('$\\sigma^2$ = ');
legend(strcat(str,num2str((sig_v.^2).','%.e')),'NumColumns',2,'Interpreter','latex');
ylabel("$10\log{\bf (E ||x_\textit{f}-x||^2/||x_\textit{f}||^2)}$",'Interpreter','latex');
xlabel('$\bf ||x_0||_2$','Interpreter','latex');
set(gca,'FontSize',16,'FontWeight','normal');
%str = sprintf('OMP n=%d, m=%d, p=%d, NTr = %d, S=%d',n,m,p,NTr,s);
%title(str,'FontSize',16,'FontWeight','normal');
xlim([0, fct(end)]); ylim([lwbnd(lns)-0.8, max(NMSE1(1,:,1))+0.8]);
% Add New Legend
a2 = axes('position',get(gca,'position'),'visible','off'); ylim([1, 10]); hold on;
p1 = plot(fct, 10*ones(numel(fct),1)*lwbnd(4),'-k','DisplayName','$s$=4','LineWidth',2);
p2 = plot(fct, 10*ones(numel(fct),1)*lwbnd(4),'-.k','DisplayName','$s$=10','LineWidth',2);
p3 = plot(fct, 10*ones(numel(fct),1)*lwbnd(4),'--k','DisplayName','Lower Bound','LineWidth',2);
legend(a2,[p1,p2,p3],'Interpreter','latex','FontSize',16,'NumColumns',3);
% use exportgraphics(gcf,'filename.pdf','ContentType','vector'); to save as
% pdf when there are multiple axes