% Original Random Sampling Algorithm [Average Sparsity]
function [Sopt, Wopt, Fopt, Fwopt] = RandSamp_Aopt_1(R,B,t,d)
    % A, B ; LDS Matrices
    % t : # time steps
    % d : sparsity level
    [n,m] = size(B);
    e_t = 1e-8;
    % R = CtrlMatrix(A,B,t);
    Sopt = [];
    Wopt = [];
    Fopt = Inf; Fwopt = Inf;
    
    for itr = 1:10
        sig = zeros(t,m); % Actuator Schedule at a time step (bool matrix)
        w = zeros(t,m); % Weighting for each actuator
        S = []; % Support Set (represent column indices)
        V = 1:t*m; % Valid Actuators
        L = pinv(R*R.');
        Prob = zeros(1,t*m); % Probabilities for actuators
        for i=1:t*m
            Prob(i) = trace(L*R(:,i)*(R(:,i).'))/n;
        end
        Prob = Prob/sum(Prob);
        CProb = cumsum([0, Prob(:).']);
        
        i = 0;
        while i < d*t
            [l,~] = histcounts(rand,CProb);
            p = V(logical(l));
            k = ceil(p/m); j = p-(k-1)*m; % Find the time step and actuator no.
            sig(k,j)=1;
            S = union(S,p);
            w(k,j) = sqrt(w(k,j).^2 + 1/(d*m*Prob(logical(l))));
            i = i+1;
        end
        w = w.';
        W = w(w>0);
        % W = reshape(w.',[],1);
        W = sqrt(t*d)*W/sqrt(sum(W.^2));
        W_Sw = R(:,S)*(diag(W).^2)*R(:,S).';
        W_S = R(:,S)*R(:,S).';
        
        fval = trace(inv(W_S + e_t*eye(n))); % A-Optimality Criteria
        if fval < Fopt
            Sopt = S;
            Wopt = W;
            Fopt = fval;
        end
        
        fwval = trace(inv(W_Sw + e_t*eye(n)));
        if fwval < Fwopt
            Fwopt = fwval;
        end
    end
end

function [V, P] = UpdateVP(V, P, sig, k)
    [t,m] = size(sig);
    for i = 1:m
        if sig(k,i) == 0
            p = (k-1)*t + i;
            idx = V~=p;
            V = V(idx); P = P(idx);
        end
    end
    P = P./sum(P);
end
