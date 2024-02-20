% Best best d actuators in each time step based on the leverage score /
% importance. Importance Sampling
function [Sopt, Fopt] = ImpSamp_Aopt(R,m,t,d)
    % A, B ; LDS Matrices
    % t : # time steps
    % d : sparsity level
    n = size(R,1);
    e_t = 0.001;
    % R = CtrlMatrix(A,B,t);
    Sopt = [];
    % Wopt = [];
    Fopt = Inf;
    
    sig = zeros(t,m); % Actuator Schedule at a time step (bool matrix)
    w = zeros(t,m); % Weighting for each actuator
    S = []; % Support Set (represent column indices)
    V = 1:t*m; % Valid Actuators
    L = pinv(R*R.');
    Prob = zeros(1,t*m); % Probabilities for actuators
    for i=1:t*m
        Prob(i) = trace(L*R(:,i)*(R(:,i).'))/n; % Importance Probability
    end
    Prob = Prob/sum(Prob);
    
    for i=0:t-1
        IProb = Prob(i*m+1:(i+1)*m); 
        [~,SAct] = sort(IProb,'descend'); % Sort the actuators acc. leverage score
        BAct = i*m + SAct(1:d); % Take best d actuators
        Sopt = union(Sopt,BAct);
    end
    
    Fopt = trace(inv(R(:,Sopt)*R(:,Sopt).'));
end