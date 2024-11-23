%% Result 9 - A-Opt vs # SubNets (Scheduler Comparision) Random Network Cluster
clear; clc; close all
rng(0);
SubNets = 1:10; NumSubNets = numel(SubNets);
s = 1; NSys = 10; c=1;
Giopt = zeros(NSys,NumSubNets);
Rwiopt = zeros(NSys,NumSubNets);
Riopt = zeros(NSys,NumSubNets);
Swiopt = zeros(NSys,NumSubNets);
Siopt = zeros(NSys,NumSubNets);
tic;
for k = 1:NumSubNets
    p = SubNets(k);
    MA = RandomCluster(p,NSys); B = eye(2*p); m =2*p;
    t = c*2*p;
    for i=1:NSys
        A = MA(:,:,i); R = CtrlMatrix(A,B,t);
        e_0 = 1e-7;
        [IW_S,~,S_ki] = FullBLI(A,B,t,s,e_0);
        [S_g,Giopt(i,k)] = GreedyScheduling_Aopt_FullB(R,IW_S,m,S_ki,t,s);
        e_0 = 1e-30;
        [S_rw,~,~,Rwiopt(i,k)] = RandSamp_Aopt(R,m,t,s,e_0);
        [S_r,Riopt(i,k)] = RandSamp_Aopt_2(R,m,t,s,e_0);
    end
end

Gopt = sum(Giopt,1)/NSys; Ropt(1,:) = sum(Riopt,1)/NSys; Ropt(2,:) = sum(Rwiopt,1)/NSys;
toc;
%% Plotting - Main Result
figure();
semilogy(SubNets,Ropt(2,:),'r--','LineWidth',3,'Marker','d','MarkerSize',10,'DisplayName','Random (Weighted)');
grid on; hold on
%semilogy(S,Sopt(3,:),'LineWidth',3,'Marker','+','MarkerSize',10,'DisplayName','Deterministic (Wo Replacement)','Color',"#000099");
semilogy(SubNets,Ropt(1,:),'r-','LineWidth',3,'Marker','d','MarkerSize',10,'DisplayName','Random (Unweighted)');
%semilogy(S,Sopt(1,:),'LineWidth',3,'Marker','+','MarkerSize',10,'DisplayName','Deterministic (Unweighted)','Color',"#009900");
semilogy(SubNets,Gopt(1,:),'b-','LineWidth',3,'Marker','s','MarkerSize',10,'DisplayName','Greedy');
%semilogy(S,Sopt(2,:),'--','LineWidth',3,'Marker','+','MarkerSize',10,'DisplayName','Deterministic (Weighted)','Color',"#009900");
%semilogy(S,Gopt(2,:),'LineWidth',3,'Marker','s','MarkerSize',10,'DisplayName','Reverse Greedy','Color',"#7E2F8E");
%semilogy(S,Lthrsh*ones(1,lg),'k-.','LineWidth',3,'DisplayName','No Sparisty Constriant');
set(gca,'FontSize',20,'FontWeight','bold')
legend();
xlabel('Number of Subnetworks $\rm{(p)}$','Interpreter','latex','FontWeight','bold','FontSize',20);
ylabel('$\rm{Tr({W_S}^{-1})}$','Interpreter','latex','FontWeight','bold','FontSize',20);
% title(['NTrails = ', num2str(NSys), ' N = ',num2str(n), ' M = ',num2str(m)])
str = sprintf('$s=%d, c=%d$, NSys=%d',s,c,NSys);
title(str,'Interpreter','latex');
%ylim([Lthrsh*0.9 max([Gopt(1), Ropt(1,1), Ropt(2,1),Sopt(2,1), Sopt(1,1)])])
xlim([SubNets(1) SubNets(end)]);
lim = max(Ropt(1,:));
ylim([1,lim*10]);