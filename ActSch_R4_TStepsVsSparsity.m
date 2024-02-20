%% Result : Greedy Algorithm Performance - TimeSteps vs Sparsity
clear; clc; close all
n = 100; m = 100;
NSys = 1;

S = 3:m/20:m/4;
lg = length(S);

%{
p = 2*log10(n)/n;
MskUt = logical(triu(ones(n),1)); % Upper Traingle Mask
Slt = binornd(1,p,n*(n-1)/2,NSys); % Bernoulli Distributed Random Numbers
MW = zeros(n,n,NSys); % Adjacency Matrices
MD = zeros(n,n,NSys); % Degree Matrices
for i=1:NSys
    Msk = zeros(n); Wi = zeros(n);
    Msk(MskUt) = logical(Slt(:,i));
    %Wi(logical(Msk)) = rand(sum(Msk,'all'),1);
    Wi(logical(Msk)) = 1;
    MW(:,:,i) = Wi+Wi.';
    MD(:,:,i) = diag(sum(Msk+Msk.'));
end
%}
MA = Erdos_Renyi(n,NSys);
%MB = rand(n,m,NSys);
MB = eye(n);
I = eye(n);

TSteps = zeros(1,lg);
TAnlys = zeros(1,lg); % row 1 - Approx bound 
TAnlys2 = zeros(1,lg); % row 2 - Lower Bound
tic;
% A = I - (MD(:,:,1)-MW(:,:,1))/n;
A = MA(:,:,1);
B = MB(:,:,1);
R = CtrlMatrix(A,B,n);
e_0 = eigs(R*R.',1,'smallestabs')/(n-1);
%e_0 = 1;
parfor k=1:lg
    s=S(k);
    [~,TSteps(k),~,TAnlys(k),TAnlys2(k)] = GreedyScheduling_Aopt_1(R,B,n,s,e_0); 
end
toc;

%% Plotting
figure();
plot(S,TAnlys,'--sr','DisplayName','Approx Bound','LineWidth',3,'MarkerSize',10)
hold on; grid on
plot(S,TSteps,'-ok','DisplayName','Empirical','LineWidth',3,'MarkerSize',12)
plot(S,TAnlys2,'-.db','DisplayName','Lower Bound','LineWidth',3,'MarkerSize',10)
h1 = legend();
set(h1,'Interpreter','latex','FontSize',36)
set(gca,'FontSize',20,'FontWeight','bold')
xlabel('Sparsity (s)')
ylabel('Time Steps')
tle = sprintf('\\epsilon_0 = %.2f N = %d M = %d',e_0,n,m);
title(tle)