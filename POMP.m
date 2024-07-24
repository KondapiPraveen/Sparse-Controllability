% POMP Sparse Recovery Algorithm
%{
    Model: b = D*sk + e     sk: Sparse Input
    D: Measurement Matrix   b: Measurements
    K: # Sparse Vectors   s: Piecewise Sparsity
%}
function sk = POMP(D,b,K,s)
    % Initialize
    tol = eps; m = size(D,2)/K; % Length of each Sparse Vector
    % Support Set ; Abandon Set
    SuppS = []; % SuppA = [];
    SetU = 1:K*m; SetV = SetU; r = b;
    S = zeros(K,m); % Binary Matrix - row : time index, col : act index
    sk = zeros(K*m,1); % solution
    Dnorm = vecnorm(D);
    %D = D./Dnorm;
    % Iteration
    for k=1:K*s
        % SetV = setdiff(SetU,SuppA);
        Dk = D(:,SetV);
        hk = Dk.'*r; % residual estimate
        [~,ij] = max(abs(hk./(Dnorm(SetV).')));
        j = SetV(ij);
        SuppS = union(SuppS, j);
        bck = ceil(j/m); p = mod(j,m)+1; % time index, act index
        S(bck,p) = 1;
        SetV = setdiff(SetV,j);
        if sum(S(bck,:)) == s
            Sbck = (bck-1)*m+1:bck*m;
            % Srst = setdiff(Sbck,intersect(SuppS,Sbck));
            % SuppA = union(SuppA,Sbck);
            SetV = setdiff(SetV,Sbck);
        end
        sk(SuppS) = pinv(D(:,SuppS))*b;
        r = b-D*sk;
        
        %{
        if norm(r)^2 < tol || numel(SuppS) == K*s
            break
        end
        %}
    end
end