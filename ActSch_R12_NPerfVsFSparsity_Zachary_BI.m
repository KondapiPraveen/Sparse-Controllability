%% Result 4 : Normalized Tr(W_S^{-1}) vs (s/m) for all actuator schedulers
clear; close all; clc; rng(0); % For Reproducability
FS = 0.1:0.1:1;
FS = union(FS,1); lg = length(FS);

e_0 = 1e-6; e_01 = 1e-20;
load ../Data/zachary.mat
zachary_lap = diag(sum(zachary_bin,2)) - zachary_bin;
n = size(zachary_bin,1); m=n; A = eye(n) - zachary_lap/n; B = eye(m);
S = floor(FS*m); t=ceil(n/3); % time to control
n_algs = 5; % Number of actuator schedulers
% 1 - RBn-Greedy, 2 -Deterministic weighted, 3 - Deterministic unweighted,
% 4 - Random unweighted, 5 - Random weighted.
Aopt = zeros(n_algs, lg);
Refopt = zeros(1,lg);

tic;
R = CtrlMatrix(A,B,t);
for k=1:lg
    s = S(k);
    Refopt(k) = trace(inv(R*R.'));
    if s == m
        Aopt(1,k) = Refopt(k);
    else
        [IW_S,~,S_ki] = FullBLI(A,B,t,s,e_0);
        [~,Aopt(1,k)] = GreedyScheduling_Aopt_FullB(R,IW_S,m,S_ki,t,s);
    end
    [S_s,~,Aopt(3,k),Aopt(2,k),~] = SparseScheduling(R,m,t,s,e_01);
    [S_r,Aopt(4,k)] = RandSamp_Aopt_2(R,m,t,s,e_01);
    %[S_rw,~,~,Aopt(5,k)] = RandSamp_Aopt(R,m,t,s,e_01);
end
toc;
Aopt = Aopt./Refopt;
%% 
load('./sparse-control/exp/R12_Zachary_BI.mat')
%s_greedy_cost = [218.5456, 93.0436, 39.0217, 25.0047, 19.9875, 17.5691, 13.9622, 12.5393, 9.1032, 8.4112];
%s_greedy_mcmc_cost = [95.1536, 33.8757, 18.7669, 14.7879, 11.9836, 10.7824, 9.9526, 9.1489, 8.7187, 8.4112];
s_greedy_cost = s_greedy_cost./Refopt; 
s_greedy_mcmc_cost = s_greedy_mcmc_cost./Refopt;
%% Plotting
figure();
x_axis = S/m;
%loglog(x_axis,Aopt(5,:),'r-','LineWidth',2,'Marker','d','MarkerSize',10,'DisplayName','Random (Weighted)')
%loglog(x_axis,Aopt(2,:),'LineWidth',2,'Marker','+','MarkerSize',10,'DisplayName','Deterministic (Weighted)','Color',"#009900")

loglog(x_axis,Aopt(3,:),'--','LineWidth',2,'Marker','+','MarkerSize',10,'DisplayName','Deterministic','Color',"#009900")
grid on; hold on;
loglog(x_axis,Aopt(4,:),'r--','LineWidth',2,'Marker','d','MarkerSize',10,'DisplayName','Random')
loglog(x_axis,s_greedy_cost,'LineWidth',2,'Marker','*','MarkerSize',10,'DisplayName','$s$-sparse greedy','Color','#A2142F')
loglog(x_axis,s_greedy_mcmc_cost,'LineWidth',2,'Marker','x','MarkerSize',10,'DisplayName','$s$-sparse greedy + mcmc','Color','#7E2F8E')
loglog(x_axis,Aopt(1,:),'b-','LineWidth',2,'Marker','s','MarkerSize',10,'DisplayName','RB$n$-greedy')
loglog(x_axis,m./S,'-.k','LineWidth',2,'Marker','o','MarkerSize',10,'DisplayName','$\frac{m}{s}$')
legend('NumColumns',1,'Interpreter','latex','Location','southwest');
set(gca,'FontSize',16,'FontWeight','normal')
xlabel('Fraction of Active Actuators per Time ($\frac{s}{m}$)','Interpreter','latex','FontWeight','normal','FontSize',17);
ylabel('$\bf E \rm {\frac{Tr({W_S}^{-1})}{Tr({W}^{-1})}}$','Interpreter','latex','FontWeight','normal','FontSize',17);
xlim([FS(1)-0.02 1.05]);
ylim([0.9, 12])
xticks([0.1, 0.2, 0.3, 0.5, 0.8, 1])
yticks([0.7, 1, 2, 3, 5, 10])