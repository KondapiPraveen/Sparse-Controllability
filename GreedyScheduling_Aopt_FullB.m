% New and Modified Greedy Actuator Scheduling
% Inputs : R - Controllabilty Matrix, m - input vector dimension
%          S - Initial Schedule Binary Lookup Table
%          ts - Control Horizon, s - sparsity level
% Outputs : S - Actutator Scheudle, Fopt - Metric Value, t - # Iterations
%           Talys - Lower Bound on t, Talys2 - Approx Bound on t
function [S, Fopt, LBnd, UBnd, Talys, Talys2] = GreedyScheduling_Aopt_FullB(R,m,S_ki,ts,s)
    n = size(R,1);
    
    % Controllability matrix and Its gramiam
    % R = CtrlMatrix(A,B,ts);
    % W_S = R(:,S)*R(:,S).';
    % IW_S = (1/e_t)*eye(n);
    % IW_S = inv(R(:,S)*R(:,S).' + e_t*eye(n));
    
	% Greedy Algorithm
    S_k = zeros(ts,m); % Track support at each time instant (Boolean matrix)
    K = size(S_ki,1);
    S_k(1:K,:) = S_ki;
    Init_Supp = m*diag(0:K-1)*S_ki+S_ki*diag(1:m);
    % 1 means the actutator is active in that time step (given by row no.)
    % Column no. implies the actuator index
    S = Init_Supp(S_ki>0); % Total Support
    V = 1:m*ts; % Valid Support
    V = setdiff(V,S);
    r=0;
    e_0 = 1e-20;
    IW_S = inv(R(:,S)*R(:,S).'+e_0*eye(n));
    while ~isempty(V)
        % optimize for the schedule
	    % p is index obtained after optimization
        % p = 0; % To track optimal actuator at a given iteration
	    % Tropt = trace(IW_S); % Optimal Trace Value
            
        v = R(:,V);
	    X = v.'*IW_S;
        Y = sum(X.*(v.'),2);
	    % [Trc, l] = max(diag(X*X.')./(1+diag(Y)));
        UpVec = (vecnorm(X,2,2).^2)./(1+Y); % Update Vector
        %{
	    [Trc, l] = max(UpVec(end:-1:1));
        l = length(V)-l+1;
        %}
	    [Trc, l] = max(UpVec);
        % [Trc, l] = max(diag(v.'*IW_S*IW_S*v)./(1+Y));
        p = V(l);
        if Trc <= 0
            fprintf('Error \n');
        end
        %{
	    for l=V
            % S_i = union(S,l); % Intermediate Support
            v = R(:,l);
            % W_Si = W_S + v*v.';
            IW_Si = IW_S - (IW_S*(v*v.')*IW_S)/(1+v.'*IW_S*v);
            % IW_Si = inv(R(:,S)*R(:,S).' + v*v.' + e_t*eye(n));
            %{
	        if norm(IW_Si - IW_Si2) > 1
                fprintf('Alert 1 \n');
            end
            %}
	        Trc = trace(IW_Si);
            %{
            if Trc < 0
    	        fprintf('Alert 2 \n');
            end
            %}
            if Trc <= Tropt
                Tropt = Trc;
                p = l;
            end
        end
        %}
        % Valid Support Update
	    V = setdiff(V,p); % Remove the selected actuator
        k = ceil(p/m); j = p-(k-1)*m; % Find the time step and actuator no.
        if sum(S_k(k,:)) < s % Check if the new actutator can be added
	        S_k(k,j)=1;
            S = union(S,p);
            v1 = R(:,p);
            IW_S = IW_S - ((IW_S*(v1*v1.')*IW_S)/(1+v1.'*IW_S*v1));
	        % IW_S = inv(R(:,S)*R(:,S).' + e_t*eye(n));
        end
        if sum(S_k(k,:)) == s % Remove other actuators in the same time step
            forbid = (k-1)*m+(1:m);
            V = setdiff(V,forbid);
        end
        r = r+1;
            
        %{
        W_S = R(:,S)*R(:,S).' + e_t*eye(n);
        IW_SI = inv(W_S);
        %}
    end
    W_S = R(:,S)*R(:,S).';
    %IW_S = inv(W_S); 
    e_t=0;
    [~,D,~] = svd(W_S + e_t*eye(n));
    Fopt = sum(1./diag(D));
    Talys = 0; Talys2 = 0;
    %{
    if t==0
        Bnd = 2*trace(IW_S) + (n)/(eigs(R*R.',1) + e_t);
        Talys = max(ceil(log(Bnd*e_t)/log(c)) + 1,1);
        Talys2 = max(ceil(log(trace(inv(R*R.'+e_t*eye(n)))*e_t)/log(c)),1);
    end
    %}
    %trace(inv(W_S));
    %Fopt = trace(IW_S);
    %fprintf('------------------------------\n')
    % fprintf('The Rank of the Sparse Scheduled Controllability Gramian is %d (Greedy Scheduling A Optimality) \n',rank(W_S))
    % fprintf('The Condition no. of the Controllability Gramian is: %d (Greedy Scheduling A Optimality) \n',cond(W_S))
    
    % Bounds
    LBnd=0; UBnd=0;
    %{
    lmu_1 = trace(W_S); lmu_2 = norm(W_S,'fro')^2; lbt = eigs(W_S,1); laph = eigs(W_S,1,'smallestabs');  
    LBnd1 = ((lmu_1*lbt - n*(lbt^2))*n + lmu_2*n - lmu_1^2)/(lbt*(lmu_2-lmu_1*lbt));
    UBnd2 = ((lmu_1*laph - n*(laph^2))*n + lmu_2*n - lmu_1^2)/(laph*(lmu_2-lmu_1*laph));
    W = R*R.';
    umu_1 = trace(W); umu_2 = norm(W,'fro')^2; ubt = eigs(W,1); uaph = eigs(W,1,'smallestabs');
    UBnd1 = ((umu_1*uaph - n*(uaph^2))*n + umu_2*n - umu_1^2)/(uaph*(umu_2-umu_1*uaph));
    LBnd2 = ((umu_1*ubt - n*(ubt^2))*n + umu_2*n - umu_1^2)/(ubt*(umu_2-umu_1*ubt));
    LBnd = LBnd1/UBnd1;
    UBnd = UBnd2/LBnd2;
    %}
end
