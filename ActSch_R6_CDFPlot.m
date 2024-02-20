%% Scheduling Result 6 : CDF Plot for Time Static and Time Varying Support
clear; clc; close all
N = 20; m = 20; % State and Input Dimensions
rng(0) % For Reproducability
NSys = 500; % # Independent Trials
e_0 = 1e-2;
%e_0 = eps;

S = [5,10,15,m];
lg = length(S);
Giopt = zeros(lg,NSys); % Time Varying Support
GSiopt = zeros(lg-1,NSys); % Fixed Support
tic;

n = N; m = n;
MA = Erdos_Renyi(n,NSys);
B = eye(m); % m=n
    
parfor l=1:NSys
    A = MA(:,:,l); %B = MB(:,:,l);
    R = CtrlMatrix(A,B,n); % Controllabilty Matrix
    for k=1:lg
        s=S(k);
        if s==m
            Giopt(k,l) = trace(inv(R*R.'));
        else
            [~,~,Giopt(k,l)] = GreedyScheduling_Aopt_1(R,m,n,s,e_0);
            if s>5
                [Supp,~,GSiopt(k,l)] = GreedyScheduling_Static_Aopt_1(R,m,n,s,e_0);
            end
        end
    end
end
toc;
%% Plotting
h = cdfplot(10*log10(Giopt(1,:))); % TV s=5
set(h,'color','r','LineStyle','-','LineWidth',3)
hold on; grid on;
h = cdfplot(10*log10(Giopt(2,:))); % TV s=10
set(h,'color','g','LineStyle','-','LineWidth',3)
h = cdfplot(10*log10(Giopt(3,:))); % TV s=15
set(h,'color','b','LineStyle','-','LineWidth',3)
h = cdfplot(10*log10(Giopt(4,:))); % TV s=200
set(h,'color','k','LineStyle','-','LineWidth',3)
h = cdfplot(10*log10(GSiopt(2,:))); % F s=10
set(h,'color','g','LineStyle','--','LineWidth',3)
h = cdfplot(10*log10(GSiopt(3,:))); % F s=15
set(h,'color','b','LineStyle','--','LineWidth',3)
xlabel('$10\log(Tr(W_S^{-1}))$','Interpreter','latex')
ylabel('CDF')

dstr = 'Dynamic Support s = '; str = strcat(dstr,num2str(S(1)));
h1 = plot(median(10*log10(Giopt(1,:))),0.5,'r-*','DisplayName',str,'LineWidth',2,'MarkerSize',8);
str = strcat(dstr,num2str(S(2)));
h2 = plot(median(10*log10(Giopt(2,:))),0.5,'g-s','DisplayName',str,'LineWidth',2,'MarkerSize',8);
str = strcat(dstr,num2str(S(3)));
h3 = plot(median(10*log10(Giopt(3,:))),0.5,'b-d','DisplayName',str,'LineWidth',2,'MarkerSize',8);
str = strcat(dstr,num2str(S(4)));
h4 = plot(median(10*log10(Giopt(4,:))),0.5,'k-o','DisplayName',str,'LineWidth',2,'MarkerSize',8);
fstr = 'fixed Support s = '; str = strcat(fstr,num2str(S(2)));
h5 = plot(median(10*log10(GSiopt(2,:))),0.5,'g--s','DisplayName',str,'LineWidth',2,'MarkerSize',8);
str = strcat(fstr,num2str(S(3)));
h6 = plot(median(10*log10(GSiopt(3,:))),0.5,'b--d','DisplayName',str,'LineWidth',2,'MarkerSize',8);
legend([h4,h3,h2,h5,h1,h6]);