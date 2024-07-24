function [S, c, Fopt1, Fwopt, Fopt2] = SparseScheduling(R,m,t,s)
    % Results from this are labeled as Deterministic Scheduling
    % A,B are System Parameters
    % t : # Time Steps
    % s : Sparsity Level at each Time Step
    % c = 0; Fwopt = 0;
    n = size(R,1);
    e_t = 1e-20;
    % R = CtrlMatrix(A,B,t);
    V = ((R*R.')^(-1/2))*R;
    U = V; % Alg 2
    % U = eye(m*t); % Alg 3
    % I = eye(m); U = (1/sqrt(t))*repmat(I,[1 t]); % Alg 4
    
    [c, S] = DualSet(V,U,s,t);
    %c = sqrt(c); % S_i(k) % Alg 3, 4
    c = sqrt(c/(1+(n/(s*t)))); % Alg 2
    c = sqrt(s*t)*c/sqrt(sum(c.^2)); % Normalization
    W_Sw = R*(diag(c)^2)*R.';
    %Fwopt = trace(inv(W_Sw + e_t*eye(n)));
    Dw = svd(W_Sw);
	Rnk = rank(W_Sw);
    if Rnk<n
        Dw(Rnk+1:end) = e_t;
    end
    Fwopt = sum(1./(Dw + e_t));
    W_S = R(:,S)*R(:,S).';
    %Fopt1 = trace(inv(W_S + e_t*eye(n)));
    D = svd(W_S);
	Rnk = rank(W_S);
    if Rnk<n
        D(Rnk+1:end) = e_t;
    end
    Fopt1 = sum(1./(D + e_t));
    
    % UnWeighted Scheduling (actuators without replacement)
    Fopt2 = 0;
    %{
    S2 = DualSet_2(V,U,s,t);
    W_S = R(:,S2)*R(:,S2).';
    Fopt2 = trace(inv(W_S + e_t*eye(n)));
    %}
end
