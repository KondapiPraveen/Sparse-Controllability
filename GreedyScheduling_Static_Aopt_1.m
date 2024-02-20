% New and Modified Greedy Actuator Scheduling
% Inputs : R - Controllabilty Matrix, m - input vector dimension
%          ts - Control Horizon, s - sparsity level
% Outputs : S - Actutator Scheudle, Fopt - Metric Value, t - # Iterations
%           Talys - Lower Bound on t, Talys2 - Approx Bound on t
function [S, t, Fopt,Talys,Talys2] = GreedyScheduling_Static_Aopt_1(R,m,ts,s,e_0)
    n = size(R,1);
    c = 2; e_t = e_0;
    
    % Controllability matrix and Its gramiam
    % R = CtrlMatrix(A,B,ts);
    t = 0;
    % W_S = R(:,S)*R(:,S).';
    % IW_S = (1/e_t)*eye(n);
    % IW_S = inv(R(:,S)*R(:,S).' + e_t*eye(n));
    CanItr = 1;
    X = Fixed_Actuators(R,ts,m);
    
    while CanItr
        % Greedy Algorithm
        V = 1:m; % Valid Support
        % 1 means the actutator is active in that time step (given by row no.)
        % Column no. implies the actuator index
        S = []; % Total Support
        r=0;
        W_S = e_t*eye(n);
        while ~isempty(V)
            % optimize for the schedule
            % p is index obtained after optimization
            % p = 0; % To track optimal actuator at a given iteration
            % Tropt = trace(IW_S); % Optimal Trace Value
            
            NValid = numel(V);
            TrVec = zeros(1,NValid);
            for act_idx = 1:NValid
                TrVec(act_idx) = trace(inv(W_S+X(:,:,V(act_idx))));
            end
            [Trc,j] = min(TrVec);
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
            p = V(j);
            V = setdiff(V,p); % Remove the selected actuator
            if numel(S) < s % Check if the new actutator can be added
                S = union(S,p);
                W_S = W_S + X(:,:,p);
                % IW_S = inv(R(:,S)*R(:,S).' + e_t*eye(n));
            end
            if numel(S) == s % Remove other actuators in the same time step
                break;
            end
            r = r+1;
            
            %{
            W_S = R(:,S)*R(:,S).' + e_t*eye(n);
            IW_SI = inv(W_S);
            %}
        end
        % Compute IW_S
        SAct = (0:ts-1)*m + S.'; % Selected Actuators
        Supp = SAct(:);
        W_S = R(:,Supp)*R(:,Supp).';
        Trcopt = trace(inv(W_S+e_t*eye(n)));
        Talys = 0; Talys2 = 0;
        %{
        if t==0
            Bnd = 2*trace(IW_S) + (n)/(eigs(R*R.',1) + e_t);
            Talys = max(ceil(log(Bnd*e_t)/log(c)) + 1,1);
            Talys2 = max(ceil(log(trace(inv(R*R.'+e_t*eye(n)))*e_t)/log(c)),1);
        end
        %}
        CanItr = Trcopt > 1/e_t;
        t = t+1;
        e_t = e_t/c;
        %trace(inv(W_S));
    end
    e_a = 0;
    Fopt = trace(inv(W_S + e_a*eye(n)));
    %fprintf('------------------------------\n')
    % fprintf('The Rank of the Sparse Scheduled Controllability Gramian is %d (Greedy Scheduling A Optimality) \n',rank(W_S))
    % fprintf('The Condition no. of the Controllability Gramian is: %d (Greedy Scheduling A Optimality) \n',cond(W_S))
end


% Input : R - Controllability Matrix, ts : time steps, m : input dimension
% Output : X(:,:,i) = [A^{k-1}B_i A^{k-2}B_i ... B_i]*(.)^T
function X = Fixed_Actuators(R,ts,m)
    n = size(R,1);
    X = zeros(n,n,m);
    supp = (0:ts-1)*m;
    for i=1:m
        Y = R(:,supp+i);
        X(:,:,i) = Y*Y.';
    end
end