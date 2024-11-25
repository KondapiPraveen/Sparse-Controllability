%% Result I : f_opt vs Sparsity - Comparision of Greedy vs Actuator Schedulers
clear; clc; close all
n = 20; m = n; % State, Input dimension, Control Time Steps
%lowlvl = 3;
rng(0);  NSys = 1; %stp = 2;

cntExample = 1; % For Counter Example
if cntExample
    load ./sparse-control/Ipexp/CounterExample.mat
    MA = A; B= 10*B;
else
   NSys = 100;
   MA = Erdos_Renyi(n,NSys);
   MB = rand(n,m,NSys);
end

lowlvl = max(n-rank(MA),2);
S = lowlvl:lowlvl+6; % Sparsity Level
lg = length(S);
Ropt = zeros(2,lg); % Row 1 - Unweighted, Row 2 - Weighted
Sopt = zeros(3,lg); % Row 1 - Unweighted, Row 2 - Weighted
Giopt = zeros(NSys,lg); % Forward Greedy
Giopt2 = zeros(NSys,lg); % Reverse Greedy
Riopt = zeros(NSys,lg); % Unweighted
Rwiopt = zeros(NSys,lg); % Weighted
Siopt2 = 1e30*ones(NSys,lg); % Unweighted (Without Replacement)
Siopt = 1e30*ones(NSys,lg); % Unweighted (From Weighted)
Swiopt = 1e30*ones(NSys,lg); % Weighted

% MA = randn(n,n,NTrails); MB = randn(n,m,NTrails);

%{
p = 2*log10(n)/n;
MskUt = logical(triu(ones(n),1)); % Upper Traingle Mask
Slt = binornd(1,p,n*(n-1)/2,NSys); % Bernoulli Distributed Random Numbers
MW = zeros(n,n,NSys); % Adjacency Matrices
MD = zeros(n,n,NSys); % Degree Matrices
for i=1:NSys
    Msk = zeros(n); Wi = zeros(n);
    Msk(MskUt) = logical(Slt(:,i));
    % Wi(logical(Msk)) = randn(sum(Msk,'all'),1);
    Wi(logical(Msk)) = 1;
    MW(:,:,i) = Wi+Wi.';
    MD(:,:,i) = diag(sum(Msk+Msk.'));
end
clear Msk MskUt Wi Slt
%}

% mdl = 'NM2 ';
% I = eye(n);
Lthrsh = zeros(NSys,1);
e_01 = 1e-5; e_0 = 1e-20;
tic;
for i = 1:NSys
    A = MA(:,:,i); %B = MB(:,:,i); % Random input matrix
    % A = I - (MD(:,:,i)-MW(:,:,i))/n;
    t = ceil(n/2); R = CtrlMatrix(A,B,t);
    %NrmZ = trace(inv(R*R.')); % Normalizing Constant
    for k=1:lg
        s = S(k); %t = ceil(n/s);
        %R = CtrlMatrix(A,B,t); NrmZ = 1; %NrmZ = trace(inv(R*R.'));
        [IW_S,~,S_ki] = FullBLI(A,B,t,s,e_01);
        [S_g,Giopt(i,k)] = GreedyScheduling_Aopt_FullB(R,IW_S,m,S_ki,t,s);
        %{
        if (t*s > n)
            [S_s,~,Siopt(i,k),Swiopt(i,k),Siopt2(i,k)] = SparseScheduling(R,m,t,s,e_0);
        end
        [S_rw,~,~,Rwiopt(i,k)] = RandSamp_Aopt(R,m,t,s,e_0);
        [S_r,Riopt(i,k)] = RandSamp_Aopt_2(R,m,t,s,e_0);
        Giopt(i,k) = Giopt(i,k)/NrmZ; Riopt(i,k) = Riopt(i,k)/NrmZ; Rwiopt(i,k) = Rwiopt(i,k)/NrmZ;
        Siopt(i,k) = Siopt(i,k)/NrmZ; Siopt2(i,k) = Siopt2(i,k)/NrmZ; Swiopt(i,k) = Swiopt(i,k)/NrmZ;
        %}
    end
    Lthrsh(i) = trace(inv(R*R.' + 0.001*eye(n)));
end
Gopt(1,:) = sum(Giopt,1); Gopt(2,:) = sum(Giopt2,1); % Forward and Reverse Greedy
Sopt(1,:) = sum(Siopt,1); Sopt(3,:) = sum(Siopt2,1); % Unweighted Deterministic
Ropt(1,:) = sum(Riopt,1); % Unweighted Randomized
Sopt(2,:) = sum(Swiopt,1); Ropt(2,:) = sum(Rwiopt,1);
Lthrsh = sum(Lthrsh);
Gopt = Gopt/NSys; Sopt= Sopt/NSys;
Ropt = Ropt/NSys; Lthrsh = Lthrsh/NSys;
toc;
%% Plotting - Main Result
figure();
semilogy(S,Ropt(2,:),'r--','LineWidth',3,'Marker','d','MarkerSize',10,'DisplayName','Random (Weighted)');
grid on; hold on
%semilogy(S,Sopt(3,:),'LineWidth',3,'Marker','+','MarkerSize',10,'DisplayName','Deterministic (Wo Replacement)','Color',"#000099");
semilogy(S,Ropt(1,:),'r-','LineWidth',3,'Marker','d','MarkerSize',10,'DisplayName','Random (Unweighted)');
semilogy(S,Sopt(1,:),'LineWidth',3,'Marker','+','MarkerSize',10,'DisplayName','Deterministic (Unweighted)','Color',"#009900");
semilogy(S,Gopt(1,:),'b-','LineWidth',3,'Marker','s','MarkerSize',10,'DisplayName','Greedy');
semilogy(S,Sopt(2,:),'--','LineWidth',3,'Marker','+','MarkerSize',10,'DisplayName','Deterministic (Weighted)','Color',"#009900");
%semilogy(S,Gopt(2,:),'LineWidth',3,'Marker','s','MarkerSize',10,'DisplayName','Reverse Greedy','Color',"#7E2F8E");
%semilogy(S,Lthrsh*ones(1,lg),'k-.','LineWidth',3,'DisplayName','No Sparisty Constriant');
set(gca,'FontSize',20,'FontWeight','bold')
legend();
xlabel('$\rm{Sparsity (s)}$','Interpreter','latex','FontWeight','bold','FontSize',20);
ylabel('$\rm{Tr({W_S}^{-1})}$','Interpreter','latex','FontWeight','bold','FontSize',20);
% title(['NTrails = ', num2str(NSys), ' N = ',num2str(n), ' M = ',num2str(m)])
title('Tr({W_S}^{-1}) vs Sparsity (s)')
%ylim([Lthrsh*0.9 max([Gopt(1), Ropt(1,1), Ropt(2,1),Sopt(2,1), Sopt(1,1)])])
xlim([S(1)-1 S(end)+1]);
lim = max(Gopt(1,:));
ylim([100,lim*10]);
%%
figure();
X = categorical(S);
RCE = ones(6,length(S));
RCE(2,:) = Gopt(2,:)./Gopt(1,:);
RCE(3,:) = Ropt(1,:)./Gopt(1,:);
RCE(4,:) = Ropt(2,:)./Gopt(1,:);
RCE(5,:) = Sopt(1,:)./Gopt(1,:);
RCE(6,:) = Sopt(2,:)./Gopt(1,:);
bar(X,RCE.')
xlabel('Sparsity Level (s)');
ylabel('Relative Tr({W_S}^{-1}) w.r.t Greedy Scheduling');
str = sprintf(' N = %d M = %d NSys = %0.1e',n,m,NSys);
legend('Forward Greedy', 'Reverse Greedy', 'Random (Unweighted)','Random (Weighted)','Deterministic (Unweighted)','Deterministic (Weighted)')
title(str);
%%
%{
figure();
semilogy(S,Ropt(2,:),'r--','LineWidth',3,'Marker','d','MarkerSize',10,'DisplayName','Random Scheduling (Weighted)');
grid on; hold on
semilogy(S,Ropt(1,:),'r-','LineWidth',3,'Marker','d','MarkerSize',10,'DisplayName','Random Scheduling (Unweighted)');
semilogy(S,Sopt(1,:),'-','LineWidth',3,'Marker','+','MarkerSize',10,'DisplayName','Deterministic Scheduling (Unweighted)','Color',"#009900");
semilogy(S,Gopt,'b-','LineWidth',3,'Marker','s','MarkerSize',10,'DisplayName','Greedy Scheduling');
semilogy(S,Sopt(2,:),'--','LineWidth',3,'Marker','+','MarkerSize',10,'DisplayName','Deterministic Scheduling (Weighted)','Color',"#009900");
semilogy(S,Lthrsh*ones(1,lg),'k-.','LineWidth',3,'DisplayName','No Sparisty Constriant');
set(gca,'FontSize',20,'FontWeight','bold')
legend();
xlabel('Sparsity (\it{s})','FontWeight','bold','FontSize',20);
ylabel('Tr({W_S}^{-1})','FontWeight','bold','FontSize',20);
% title(['NTrails = ', num2str(NSys), ' N = ',num2str(n), ' M = ',num2str(m)])
title('Tr({W_S}^{-1}) vs Sparsity (s)')
ylim([Lthrsh*0.9 max([Gopt(1), Ropt(1,1), Ropt(2,1),Sopt(2,1), Sopt(1,1)])])
%%
figure();
X = categorical(S);
RCE = ones(5,length(S));
RCE(2:3,:) = Ropt./Gopt;
RCE(4:5,:) = Sopt./Gopt;
bar(X,RCE.')
xlabel('Sparsity Level (s)');
ylabel('Relative Tr({W_S}^{-1}) w.r.t Greedy Scheduling');
str = sprintf(' N = %d M = %d NSys = %0.1e',n,m,NSys);
legend('Greedy Scheduling','Random Sampling','Random Sampling (Weighted)','Deterministic Scheduling','Deterministic Scheduling (Weighted)')
title(strcat(mdl,str));
%}
