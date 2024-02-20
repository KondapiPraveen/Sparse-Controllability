%% Scheduling Result 5 : |W_S|/|W| vs Sparsity Performance Analysis of A-Opt Greedy Algorithm
clear; clc; close all
N = 10:10:40; m = 20;
rng(0)
NSys = 30;
e_0 = 1e-10;
%e_0 = eps;

S = (0.1:0.1:1)*m;
lg = length(S);
ln = length(N);
Giopt_max = zeros(ln,lg,NSys); % Max E.Val
Giopt_min = zeros(ln,lg,NSys); % Min E.Val
nGopt_max = zeros(ln,lg,NSys);
nGopt_min = zeros(ln,lg,NSys);

tic;
for j=1:ln
    n = N(j);
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

    MB = rand(n,m,NSys);
    
    I = eye(n); %B = eye(n);
    parfor l=1:NSys
        A = I - (MD(:,:,l)-MW(:,:,l))/n;
        B = MB(:,:,l);
        R = CtrlMatrix(A,B,n);
        for k=1:lg
            s=S(k);
            if k==lg
                W = R*R.';
                Giopt_max(j,k,l) = eigs(W,1);
                Giopt_min(j,k,l) = eigs(W,1,'smallestabs');
            else
               [Supp] = GreedyScheduling_Aopt_1(R,m,n,s,e_0);
               W_S = R(:,Supp)*R(:,Supp).';
               Giopt_max(j,k,l) = eigs(W_S,1);
               Giopt_min(j,k,l) = eigs(W_S,1,'smallestabs');
            end
        end
    end
end
toc;

for l=1:NSys
    nGopt_max(:,:,l) = Giopt_max(:,:,l)./Giopt_max(:,end,l);
    nGopt_min(:,:,l) = Giopt_min(:,:,l)./Giopt_min(:,end,l);
end
Gopt_max_max = max(nGopt_max,[],3);
Gopt_max_min = min(nGopt_max,[],3);
Gopt_min_max = max(nGopt_min,[],3);
Gopt_min_min = min(nGopt_min,[],3);

%% Plotting
close all; figure();
plot(S/m,Gopt_max_max.','LineWidth',3,'Marker','s','MarkerSize',10);
xlim([0.1, 1]); ylim([0.1,1])
ylabel('Fraction of Active Actuators ($\frac{s}{m}$)','Interpreter','latex');

figure();
plot(S/m,Gopt_max_min.','LineWidth',3,'Marker','s','MarkerSize',10);
xlim([0.1, 1]); ylim([0.1,1])
ylabel('Fraction of Active Actuators ($\frac{s}{m}$)','Interpreter','latex');