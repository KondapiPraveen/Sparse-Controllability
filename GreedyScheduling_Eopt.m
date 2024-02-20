function [S, t, Fopt] = GreedyScheduling_Eopt(R,B,ts,s)
    [n,m] = size(B);
    c = 2; e_t = 0.05;
    
    % Controllability matrix and Its gramiam
    % R = CtrlMatrix(A,B,ts);
    S = []; % Support Set
    t = 0;
    W_S = R(:,S)*R(:,S).';
    % eigs(W_S+e_t*eye(n),1,'smallestabs')
    
    while 1/eigs(W_S + e_t*eye(n),1,'smallestabs') >= 1/e_t
        e_t = e_t/c;
        % Greedy Algorithm
        V = 1:m*ts; % Valid Support
        S_k = zeros(n,m); % Track support at each time instant (Boolean matrix)
        % 1 means the actutator is active in that time step (given by row no.)
        % Column no. implies the actuator index
        S = []; % Total Support
        r=0;
        while ~isempty(V)
            % optimize for the schedule
            % p is index obtained after optimization
            Esopt = inf; % optimal trace value;
            p = 0; % To track optimal actuator at a given iteration
            for l=V
                S_i = union(S,l); % Intermediate Support
                W_S = R(:,S_i)*R(:,S_i).';
                Es = 1/eigs(W_S + e_t*eye(n),1,'smallestabs');
                if Es < Esopt
                    Esopt = Es;
                    p = l;
                end
            end
            % Valid Support Update
            V = setdiff(V,p); % Remove the selected actuator
            k = ceil(p/m); j = mod(p,m) + 1; % Find the time step and actuator no.
            if sum(S_k(k,:)) < s % Check if the new actutator can be added
                S_k(k,j)=1;
                S = union(S,p);
            else % Remove other actuators in the same time step
                forbid = (k-1)*m+(1:m);
                V = setdiff(V,forbid);
            end
            r = r+1;
        end
        t = t+1;
        W_S = R(:,S)*R(:,S).';
    end
    Fopt = trace(inv(W_S));
    % fprintf('The Rank of the Sparse Scheduled Controllability Gramian is %d (Greedy Scheduling A Optimality) \n',rank(W_S))
    % fprintf('The Condition no. of the Controllability Gramian is: %d (Greedy Scheduling A Optimality) \n',cond(W_S))
end
