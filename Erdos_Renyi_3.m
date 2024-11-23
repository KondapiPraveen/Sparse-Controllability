% Generate Erdos-Renyi Graph
% Input : n - # nodes
% Input : p - Edge Probability
% Input : NSys - # Systems
function MA = Erdos_Renyi_3(n,p,NSys)
    MA = zeros(n,n,NSys);
    %p = 2*log(n)/n;
    MskUt = logical(triu(ones(n),1)); % Upper Traingle Mask
    Slt = binornd(1,p,n*(n-1)/2,NSys); % Bernoulli Distributed Random Numbers
    MW = zeros(n,n,NSys); % Adjacency Matrices
    MD = zeros(n,n,NSys); % Degree Matrices
    I = eye(n);
    for i=1:NSys
        Msk = zeros(n); Wi = zeros(n);
        Msk(MskUt) = logical(Slt(:,i));
        % Wi(logical(Msk)) = randn(sum(Msk,'all'),1);
        Wi(logical(Msk)) = 1;
        MW(:,:,i) = Wi+Wi.';
        MD(:,:,i) = diag(sum(Msk+Msk.'));
        MA(:,:,i) = I - (MD(:,:,i)-MW(:,:,i))/n;
        % Normalize the matrix with largest eigenvalue
        %l1 = eigs(MA(:,:,i),1); MA(:,:,i) = MA(:,:,i)/l1;
    end
    clear Msk MskUt Wi Slt
    
    ExportData = true;
    if ExportData && NSys == 1
        A = MA; B =I;
        save('./sparse-control/Ipexp/AER_BI.mat','A','B')
    end
end