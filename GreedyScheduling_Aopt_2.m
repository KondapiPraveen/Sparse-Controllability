% Reverse Greedy Actuator Scheduling
% Inputs : R - Controllabilty Matrix, m - input vector dimension
%          ts - Control Horizon, s - sparsity level
% Outputs : S - Actutator Scheudle, Fopt - Metric Value
function [S, Fopt] = GreedyScheduling_Aopt_2(R,m,ts,s)
    % Initialization
    e_0 = eps; % for stability of the algorithm
    n = size(R,1); % State vector dimension
    V = 1:m*ts; % Feasible set of actuators
    S_k = ones(ts,m); % Boolean Matrix to track the support at each time
    % 1 means the actutator is active in that time step (given by row no.)
    % Column no. implies the actuator index
    S = []; % Total Support
    W = R*R.' + e_0*eye(n); % Controllability Gramian
    IW_S = inv(W); % Inverse of the Gram
    
    % Iteration
    while ~isempty(V)
        % Finding optimal p = (i,j)
        v = R(:,V);
        X = v.'*IW_S;
        Y = sum(X.*(v.'),2);
        UpVec = (vecnorm(X,2,2).^2)./(1-Y); % Update Vector
        [~, l] = min(UpVec(end:-1:1));
        l = length(V)-l+1;
        p = V(l);
        % Update
        i = ceil(p/m); % Time Index
        j = mod(p,m) + 1; % Actuator Index
        if sum(S_k(i,:)) > s
            S_k(i,j) = 0; V = setdiff(V,p); % Remove the actutator 
            % from the support set and update the feasible set
            v1 = R(:,p); % The vector to be removed
            % Update IW_S
            IW_S = IW_S + ((IW_S*(v1*v1.')*IW_S)/(1-v1.'*IW_S*v1));
        end
        if sum(S_k(i,:)) == s
            fbd = (i-1)*m + (1:m); % Forbidden Actuators
            V = setdiff(V,fbd); % Remove the Forbidden Actuators
            % from Feasible Set V
        end
    end
    
    % Update S
    M = 1:m; % Indices of actuators
    for i=1:ts
        ActV = (i-1)*m + M(logical(S_k(i,:))); % Active Actuators at time step i
        S = union(S,ActV); % Add active actuators to Total Support
    end
    e_t = 0;
    Fopt = trace(inv(R(:,S)*R(:,S).'+e_t*eye(n)));
end