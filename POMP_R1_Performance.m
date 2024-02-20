%% Result IV : POMP Algorithm
clear; clc; close all
rng(0);
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
MA = Erdos_Renyi(n,NSys);
q = n;

% Simulation
stp = 2;
S = lowlvl:stp:15;
S = union(S,n);
del = 1e-6;
lg = length(S);
K = zeros(1,lg);
ACE = zeros(1,lg);
x0 = zeros(n,1);
% B = eye(n); NM 1
MB = randn(n,m,NSys); % NM2
I = eye(n);

tic;
for j=1:NSys
    % A = I - (MD(:,:,j)-MW(:,:,j))/n;
    A = MA(:,:,j); B = MB(:,:,j);
    R = CtrlMatrix(A,B,t);
    k=0;
    for s=S
        k=k+1;
        for i=1:NTrails
            w = 0;
            for l=del
                w=w+1;
                xf = randn(n,1); xf = xf/norm(xf); 
                x0 = randn(n,1); x0 = x0/norm(x0); % NM 2
                [u_p,Ki,~,~] = SpaIpDsg_1(x0, xf, A, B, s, q, l);
                K(w,k) = K(w,k) + Ki-1;
                ACE(k) = ACE(k) + norm(u_p)^2;
            end
        end
    end
end
K = K/NTot;
ACE = ACE/NTot;
toc;
%% Lower and Upper Bound
LBnd = ceil(n./S);
UBndi = zeros(2,lg);
UBndi(1,:) = n*ceil(n./S);
UBndi(2,:) = n-S+1;
UBnd = min(UBndi);