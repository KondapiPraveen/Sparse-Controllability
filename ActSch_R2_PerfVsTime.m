%% Result VI : Tr(Inv(W_s)) vs TimeSteps
clear; clc; close all
tic;
rng(0)
% System Parameters
n = 100; m = n/2;
NSys = 1;

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
MA = Erdos_Renyi(n,NSys);
% B = eye(n);
MB = rand(n,m,NSys);
I = eye(n);
q = n; e_0 = 0.01;

K = (1:0.5:3)*n;
% S = 4:n/10:m;
S = [25, 35, 45];
sl = length(S);
kl = length(K);
Giopt = zeros(sl, length(K));
Siopt = zeros(sl, length(K));
Riopt = zeros(sl, length(K));
Sopt = zeros(sl,length(K));
Gopt = zeros(sl, length(K));
Ropt = zeros(sl, length(K));
Lthrsh = zeros(2, length(K));

tic;
for j=1:NSys
    %A = I - (MD(:,:,j)-MW(:,:,j))/n; 
    A = MA(:,:,j); B = MB(:,:,j);
    parfor k=1:kl
        R = CtrlMatrix(A,B,K(k));
        for s=1:sl
            [S_s,~,Siopt(s,k),~] = SparseScheduling(R,B,K(k),S(s));
            [S_g,~,Giopt(s,k)] = GreedyScheduling_Aopt_1(R,B,K(k),S(s),e_0);
            [S_r,Riopt(s,k)] = RandSamp_Aopt_2(R,B,K(k),S(s));
        end
        Lthrsh(1,k) = trace(inv(R*R.' + 0.001*eye(n)));
    end
    Gopt = Gopt+Giopt;
    Sopt = Sopt+Siopt;
    Ropt = Ropt+Riopt;
    Lthrsh(2,:) = Lthrsh(2,:) + Lthrsh(1,:);
end
Gopt = Gopt/NSys;
Lthrsh(2,:) = Lthrsh(2,:)/NSys;
toc;
%% Plotting
figure();
LStyles = {'-','--','-.'};
for k=1:sl
    str = sprintf(' s = %d',S(k));
    semilogy(K,Gopt(k,:),strcat(LStyles{k},'s'),'DisplayName',strcat('Greedy ',str),'LineWidth',3,'MarkerSize',10);
    hold on; grid on
    semilogy(K,Sopt(k,:),strcat(LStyles{k},'+'),'DisplayName',strcat('Deterministic ',str),'LineWidth',3,'MarkerSize',10);
    semilogy(K,Ropt(k,:),strcat(LStyles{k},'d'),'DisplayName',strcat('Randomized ',str),'LineWidth',3,'MarkerSize',10);
end
semilogy(K,Lthrsh(2,:),'-.ko','DisplayName','No Sparsity Constraint','LineWidth',3,'MarkerSize',10);
set(gca,'FontSize',20,'FontWeight','bold')
legend();
xlabel('Time Steps (K)','FontWeight','bold','FontSize',20);
ylabel('Tr({W_S}^{-1})','FontWeight','bold','FontSize',20)
str = sprintf("NSys = %d N = %d M = %d",NSys,n,m);
title(str)