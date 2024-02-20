close all; figure();
lg = length(S);
semilogy(S,Rwopt(1,:),'r--','LineWidth',3,'Marker','d','MarkerSize',10,'DisplayName','Random (Weighted)');
grid on; hold on
semilogy(S,Ropt(1,:),'r-','LineWidth',3,'Marker','d','MarkerSize',10,'DisplayName','Random (Unweighted)');
semilogy(S,Sopt,'-','LineWidth',3,'Marker','+','MarkerSize',10,'DisplayName','Deterministic (Unweighted)','Color',"#009900");
semilogy(S,Swopt,'--','LineWidth',3,'Marker','+','MarkerSize',10,'DisplayName','Deterministic (Weighted)','Color',"#009900");
semilogy(S,Gopt,'b-','LineWidth',3,'Marker','s','MarkerSize',10,'DisplayName','Greedy Scheduling');
semilogy(S,Lthrsh*ones(1,lg),'k-.','LineWidth',3,'DisplayName','No Sparisty Constriant');
set(gca,'FontSize',20,'FontWeight','bold')
legend();
xlabel('Sparsity (\it{s})','FontWeight','bold','FontSize',20);
ylabel('Tr({W_S}^{-1})','FontWeight','bold','FontSize',20);
% title(['NTrails = ', num2str(NSys), ' N = ',num2str(n), ' M = ',num2str(m)])
title('Tr({W_S}^{-1}) vs Sparsity (s)')
ylim([Lthrsh*0.9 max([Gopt(1), Ropt(1), Rwopt(1),Swopt(1), Sopt(1)])])

axes('position',[.17,.17,.25,.25])
box on

idx = [1,2];
semilogy(S(idx),Sopt(idx),'-','LineWidth',3,'Marker','+','MarkerSize',10,'Color',"#009900");
grid on, hold on
semilogy(S(idx),Swopt(idx),'--','LineWidth',3,'Marker','+','MarkerSize',10,'Color',"#009900");
semilogy(S(idx),Gopt(idx),'b-','LineWidth',3,'Marker','s','MarkerSize',10);
set(gca,'FontSize',10,'FontWeight','bold')

axes('position',[.45,.17,.25,.2])
box on

idx = [5,6,7];
semilogy(S(idx),Sopt(idx),'-','LineWidth',3,'Marker','+','MarkerSize',10,'Color',"#009900");
grid on, hold on
semilogy(S(idx),Swopt(idx),'--','LineWidth',3,'Marker','+','MarkerSize',10,'Color',"#009900");
semilogy(S(idx),Ropt(idx),'r-','LineWidth',3,'Marker','d','MarkerSize',10);
semilogy(S(idx),Gopt(idx),'b-','LineWidth',3,'Marker','s','MarkerSize',10);
set(gca,'FontSize',10,'FontWeight','bold')