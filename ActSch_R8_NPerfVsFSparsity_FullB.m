%% Scheduling Result 3 : Normalized Tr(W_S^{-1}) vs (s/m) Performance Analysis of A-Opt Greedy Algorithm
clear; close all; clc; rng(0); % For Reproducability
% Initialization
FS = 0.1:0.1:1;
NSys = 100; % # Independent Trials

N = 60:20:100; % State, Input Dimensions % Setup - 1 

%{
% Setup - 2
N = 20:20:80; m = 100;
FS = union(FS,0.05);
%}
lg = length(FS);
ln = length(N);
Giopt = zeros(ln,lg,NSys); % Time Varying Support
LBnd = zeros(ln,lg,NSys); % 1 - Lower Bound - Validity of Representation
UBnd = zeros(ln,lg,NSys); % 2 - Upper Bound - Validity of Representation
GSiopt = zeros(ln,lg,NSys); % Fixed Support
tic;

for j=1:ln
    n = N(j); 
    m = n; % Setup - 1;
    MA = Erdos_Renyi(n,NSys);
    MB = rand(n,m,NSys);
    %B = eye(n); % m=n
    
    S = ceil(FS*m);
    parfor l=1:NSys
        A = MA(:,:,l); B = MB(:,:,l);
        R = CtrlMatrix(A,B,n);
        for k=1:lg
            s=S(k);
            if s==m
                Giopt(j,k,l) = trace(inv(R*R.'));
                %GSiopt(j,k,l) = trace(inv(R*R.'));
                W = R*R.';
                umu_1 = trace(W); umu_2 = norm(W,'fro')^2; ubt = eigs(W,1); uaph = eigs(W,1,'smallestabs');
                UBnd1 = ((umu_1*uaph - n*(uaph^2))*n + umu_2*n - umu_1^2)/(uaph*(umu_2-umu_1*uaph));
                LBnd2 = ((umu_1*ubt - n*(ubt^2))*n + umu_2*n - umu_1^2)/(ubt*(umu_2-umu_1*ubt));
                LBnd(j,k,l) = LBnd2/UBnd1;
                UBnd(j,k,l) = UBnd1/LBnd2;
            else
               [~,S_ki] = FullBLI(A,B,n,s);
               [~,Giopt(j,k,l)] = GreedyScheduling_Aopt_FullB(R,m,S_ki,n,s);
               %[~,~,Giopt(j,k,l), LBnd(j,k,l), UBnd(j,k,l)] = GreedyScheduling_Aopt_1(R,m,n,s,e_0);
               %[~,~,GSiopt(j,k,l)] = GreedyScheduling_Static_Aopt_1(R,m,n,s,e_0);
            end
        end
    end
end
toc;
Gopt = sum(Giopt,3)/NSys;
MLBnd = mean(LBnd,3);
MUBnd = mean(UBnd,3);
%GSopt = sum(GSiopt,3)/NSys;
%%
nGopt = zeros(ln,lg,NSys);
nGSopt = zeros(ln,lg,NSys);
for l=1:NSys
    nGopt(:,:,l) = Giopt(:,:,l)./Giopt(:,end,l);
    %nGSopt(:,:,l) = GSiopt(:,:,l)./GSiopt(:,end,l);
end
iGopt = sum(nGopt,3)/NSys; % average of normalized average energy across systems
%iGSopt = sum(nGSopt,3)/NSys;

iGopt2 = Gopt./Gopt(:,end); % normalized average across systems
%iGSopt2 = GSopt./GSopt(:,end);
%% Plotting - Time Varying Support
% Shut Down plotting for Time Varying Support
%
close all
figure();
semilogy(S,Gopt.','LineWidth',3,'Marker','s','MarkerSize',10);
grid on
legend(strcat('n = ',num2str(N.')),'NumColumns',2)
set(gca,'FontSize',20,'FontWeight','bold')
xlabel('Sparsity (\it{s})','FontWeight','bold','FontSize',20);
ylabel('Tr({W_S}^{-1})','FontWeight','bold','FontSize',20);
title(['Tr({W_S}^{-1}) vs Sparsity (s) m = ',num2str(m),' NTrails = ',num2str(NSys)])

figure();
loglog(S/m,iGopt.','LineWidth',3,'Marker','s','MarkerSize',10)
grid on; hold on
%plot(S/m,(m./S).*(1+exp(-S/m)),'-.k','LineWidth',3,'Marker','o','MarkerSize',10); str = "$\frac{m}{s}(1+e^{-\frac{s}{m}})$";
loglog(S/m,m./S,'-.k','LineWidth',3,'Marker','o','MarkerSize',10); str = "$\frac{m}{s}$"; 
h1 = legend([string(strcat('n = ',num2str(N.')));str],'NumColumns',2);
set(h1,'Interpreter','latex','FontSize',36)
set(gca,'FontSize',20,'FontWeight','bold')
xlabel('Fraction of Active Actuators per Time ($\frac{s}{m}$)','Interpreter','latex','FontWeight','bold','FontSize',24);
ylabel('$\bf E \rm {\frac{Tr({W_S}^{-1})}{Tr({W}^{-1})}}$','Interpreter','latex','FontWeight','bold','FontSize',30);
title(['Tr({W_S}^{-1}) vs Sparsity (s) m = ',num2str(m),' NTrails = ',num2str(NSys)])
xlim([FS(1)-0.01 1.05]);

figure();
rng(0);
DefaultColors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; ...
    0.4660 0.6740 0.1880; 0.3010 0.7350 0.9330; 0.6350 0.0780 0.1840];
Colors = repmat(DefaultColors(1:ln,:),3,1);
colororder(Colors)
loglog(S/m,iGopt2.','LineWidth',3,'Marker','s','MarkerSize',10)
grid on; hold on
%plot(S/m,(m./S).*(1+exp(-S/m)),'-.k','LineWidth',3,'Marker','o','MarkerSize',10); str = "$\frac{m}{s}(1+e^{-\frac{s}{m}})$";
loglog(S/m,MLBnd,'--','LineWidth',3,'Marker','*','MarkerSize',10) % Comment out to Remove Lower Bound from Plot
loglog(S/m,MUBnd,':','LineWidth',3,'Marker','*','MarkerSize',10) % Comment out to Remove Upper Bound from Plot
loglog(S/m,m./S,'-.k','LineWidth',3,'Marker','o','MarkerSize',10); str = "$\frac{m}{s}$";
%h1 = legend([string(strcat('n = ',num2str(N.')));str],'NumColumns',2);
%Uncomment above in case Upper and Lower Bounds are Commented Out and
%Comment the below line
h1 = legend([string(strcat('n = ',num2str(N.')));string(strcat('LBnd n = ',num2str(N.')));string(strcat('UBnd n = ',num2str(N.')));str],'NumColumns',2);
set(h1,'Interpreter','latex','FontSize',20)
set(gca,'FontSize',20,'FontWeight','bold')
xlabel('Fraction of Active Actuators per Time ($\frac{s}{m}$)','Interpreter','latex','FontWeight','bold','FontSize',24);
ylabel("$\frac{\bf E \rm \{ Tr({W_S}^{-1}) \} }{\bf E \rm \{ Tr({W}^{-1}) \} }$",'Interpreter','latex','FontWeight','bold','FontSize',30);
title(['Tr({W_S}^{-1}) vs Sparsity (s) m = ',num2str(m),' NTrails = ',num2str(NSys)])
xlim([FS(1)-0.01 1.05]);
%}

%% Plotting - Fixed Support
% Shut Down plotting for Fixed Support
%{
close all
figure();
semilogy(S,GSopt.','LineWidth',3,'Marker','s','MarkerSize',10);
grid on
legend(strcat('n = ',num2str(N.')),'NumColumns',2)
set(gca,'FontSize',20,'FontWeight','bold')
xlabel('Sparsity (\it{s})','FontWeight','bold','FontSize',20);
ylabel('Tr({W_S}^{-1})','FontWeight','bold','FontSize',20);
title(['Tr({W_S}^{-1}) vs Sparsity (s) m = ',num2str(m),' NTrails = ',num2str(NSys)])

figure();
loglog(S/m,iGSopt.','LineWidth',3,'Marker','s','MarkerSize',10)
grid on; hold on
%plot(S/m,(m./S).*(1+exp(-S/m)),'-.k','LineWidth',3,'Marker','o','MarkerSize',10); str = "$\frac{m}{s}(1+e^{-\frac{s}{m}})$";
loglog(S/m,m./S,'-.k','LineWidth',3,'Marker','o','MarkerSize',10); str = "$\frac{m}{s}$"; 
h1 = legend([string(strcat('n = ',num2str(N.')));str],'NumColumns',2);
set(h1,'Interpreter','latex','FontSize',36)
set(gca,'FontSize',20,'FontWeight','bold')
xlabel('Fraction of Active Actuators per Time ($\frac{s}{m}$)','Interpreter','latex','FontWeight','bold','FontSize',24);
ylabel('$\bf E \rm {\frac{Tr({W_S}^{-1})}{Tr({W}^{-1})}}$','Interpreter','latex','FontWeight','bold','FontSize',30);
title(['Tr({W_S}^{-1}) vs Sparsity (s) m = ',num2str(m),' NTrails = ',num2str(NSys)])
xlim([FS(1)-0.03 1.05]);

figure();
loglog(S/m,iGSopt2.','LineWidth',3,'Marker','s','MarkerSize',10)
grid on; hold on
%plot(S/m,(m./S).*(1+exp(-S/m)),'-.k','LineWidth',3,'Marker','o','MarkerSize',10); str = "$\frac{m}{s}(1+e^{-\frac{s}{m}})$";
loglog(S/m,m./S,'-.k','LineWidth',3,'Marker','o','MarkerSize',10); str = "$\frac{m}{s}$"; 
h1 = legend([string(strcat('n = ',num2str(N.')));str],'NumColumns',2);
set(h1,'Interpreter','latex','FontSize',36)
set(gca,'FontSize',20,'FontWeight','bold')
xlabel('Fraction of Active Actuators per Time ($\frac{s}{m}$)','Interpreter','latex','FontWeight','bold','FontSize',24);
ylabel("$\frac{\bf E \rm \{ Tr({W_S}^{-1}) \} }{\bf E \rm \{ Tr({W}^{-1}) \} }$",'Interpreter','latex','FontWeight','bold','FontSize',30);
title(['Tr({W_S}^{-1}) vs Sparsity (s) m = ',num2str(m),' NTrails = ',num2str(NSys)])
xlim([FS(1)-0.03 1.05]);
%}