% Unweighted Random Scheduling (Adapted to Piecewise Sparsity)
function [Sopt, Fopt] = RandSamp_Aopt_2(R,m,t,d,e_0)
    % A, B ; LDS Matrices
    % t : # time steps
    % d : sparsity level
    n = size(R,1);
    %e_t = 1e-20;
    % R = CtrlMatrix(A,B,t);
    Sopt = [];
    % Wopt = [];
    Fopt = Inf;
    
    L = pinv(R*R.');
    LvgScore = zeros(1,t*m); % Probabilities for actuators
    for i=1:t*m
        LvgScore(i) = trace(L*R(:,i)*(R(:,i).'))/n;
    end
    for itr = 1:5e1
        sig = zeros(t,m); % Actuator Schedule at a time step (bool matrix)
        w = zeros(t,m); % Weighting for each actuator
        S = []; % Support Set (represent column indices)
        V = 1:t*m; % Valid Actuators
        Prob = LvgScore;
        csum = sum(Prob);
        Prob = Prob/sum(Prob);
        CProb = cumsum([0, Prob(:).']);
        
        i = 0;
        rnval = rand((d+1)*t,1); cj = 1;
        while i < d*t
            [l,~] = histcounts(csum*rnval(cj),CProb);
            ll = logical(l);
            p = V(logical(l));
            k = ceil(p/m); j = p-(k-1)*m; % Find the time step and actuator no.
            if (sum(sig(k,:)) < d) || (sum(sig(k,:)) <= d && sig(k,j) == 1) % Check if the new actutator can be added
                sig(k,j)=1;
                S = union(S,p);
                % w(k,j) = sqrt(w(k,j).^2 + 1/(d*t*Prob(logical(l))));
                i = i+1;
                V = V(~ll); Prob = Prob(~ll);
                
            else
                [V, Prob] = UpdateVP(V,Prob,sig,k);
            end
            CProb = cumsum([0, Prob(:).']);
            csum = CProb(end);
            cj=cj+1;
        end
        % w = w.';
        % W = w(w>0);
        % W = reshape(w.',[],1);
        % W = sqrt(t*d)*W/sqrt(sum(W.^2));
        % W_Sw = R(:,S)*(diag(W).^2)*R(:,S).';
        W_S = R(:,S)*R(:,S).';
        
        %fval = trace(inv(W_S + e_t*eye(n))); % A-Optimality Criteria
        D = svd(W_S);
        Rnk = rank(W_S);
        if Rnk<n
            D(Rnk+1:end) = e_0;
        end
        fval = sum(1./(D + e_0));
        if fval < Fopt
            Sopt = S;
            % Wopt = W;
            Fopt = fval;
        end
    end
end

function [V, P] = UpdateVP(V, P, sig, k)
    [t,m] = size(sig);
    TAct = 1:m;
    p = (k-1)*m + TAct(sig(k,:)==0);
    [V, idx] = setdiff(V,p);
    P = P(idx);
    %{
    for i = 1:m
        if sig(k,i) == 0
            p = (k-1)*t + i;
            idx = V~=p;
            V = V(idx); P = P(idx);
        end
    end
    % P = P./sum(P);
    %}
end
