%% Result III : Avg_Time_Steps vs Sparsity - POMP Algorithm
clear; clc; close all
tic;
% System Parameters
n = 40; m = n; t = n;
lowlvl = 3;
fct = 1e2;
NTot = fct*1e2; NSys = fct*1;
NTrails = NTot/NSys;

% System Matrices Generation
% MA = randn(n,n,NSys); MB = randn(n,m,NSys);
%{
p = 2*log10(n)/n;
MskUt = logical(triu(ones(n),1)); % Upper Traingle Mask
Slt = binornd(1,p,n*(n-1)/2,NSys); % Bernoulli Distributed Random Numbers
MW = zeros(n,n,NSys); % Adjacency Matrices
MD = zeros(n,n,NSys); % Degree Matrices
for i=1:NSys
    Msk = zeros(n); Wi = zeros(n);
    Msk(MskUt) = logical(Slt(:,i));
    Wi(logical(Msk)) = randn(sum(Msk,'all'),1);
    MW(:,:,i) = Wi+Wi.';
    MD(:,:,i) = diag(sum(Msk+Msk.'));
end
clear Msk MskUt Wi Slt
%}
q = n;
MA = Erdos_Renyi(n,NSys);


% Simulation
tic;
stp = 2;
S = lowlvl:stp:15;
del = [1, 1e-1, 1e-2, 1e-4, 1e-8];
lg = length(S);
K = zeros(length(del),lg);
Kr = zeros(1,lg);
x0 = zeros(n,1);
B = eye(n); I = eye(n);
%{
k=0;
for s=S
    k = k+1;
    for j=1:NSys
        A = I - (MD(:,:,j)-MW(:,:,j))/n; % B = MB(:,:,j);
        % At = A^t;
        R = CtrlMatrix(A,B,t);
        for i=1:NTrails
            w = 0;
            for l=del
                w = w+1;
                xf = randn(n,1); xf = xf/norm(xf); 
                [u_p,Ki,~,~] = SpaIpDsg_1(x0, xf, A, B, s, q, l);
                K(w,k) = K(w,k) + Ki-1;
            end
        end
    end
end
%}
for j=1:NSys
    % A = I - (MD(:,:,j)-MW(:,:,j))/n;
    A = MA(:,:,j);
    R = CtrlMatrix(A,B,t);
    k=0;
    for s=S
        k=k+1;
        [S_r,W_r,~,~] = RandSamp_Aopt(R,B,t,s);
        Kr(k) = Kr(k) + ceil(max(S_r)/m);
        for i=1:NTrails
            w = 0;
            for l=del
                w=w+1;
                xf = randn(n,1); xf = xf/norm(xf); 
                [u_p,Ki,Rsd,~] = SpaIpDsg_1(x0, xf, A, B, s, q, l);
                K(w,k) = K(w,k) + Ki-1;
            end
        end
    end
end
K = K/(NTot);
Kr = Kr/NSys;
toc;
%% Lower and Upper Bound
LBnd = ceil(n./S);
UBndi = zeros(2,lg);
UBndi(1,:) = n*ceil(n./S);
UBndi(2,:) = n-S+1;
UBnd = min(UBndi);
%% Plotting
figure();
plot(S,K.','LineWidth',2.5,'Marker','o','MarkerSize',8)
legend(strcat("\delta = ",num2str(del.')));
hold on
plot(S,Kr,'--','LineWidth',3,'DisplayName','Random Scheduling','Marker','d','MarkerSize',10);
plot(S,n*ones(1,lg),'-.','LineWidth',3,'DisplayName','Greedy Scheduling','Marker','s','MarkerSize',8)
plot(S,LBnd,'--','LineWidth',2.5,'DisplayName','Lower Bound','Marker','o','MarkerSize',8)
plot(S,UBnd,'-.','LineWidth',2.5,'DisplayName','Upper Bound','Marker','o','MarkerSize',8)
xlabel('Sparsity (s)')
ylabel('Average Time Steps')
title('Average Number of Time Steps of POMP Algorithm')
ylim([0 n+2]);