% Algorithm I: Sparse Input Design based on initial and final state
function [u, K, Residue, j] = SpaIpDsg(x0, xf, A, B, s, q)
    n = length(x0); % dimension of state vector
    m = size(B,2); % dimension of input vector
    del = 1e-6;
    r = 2*del; % residue
    R_B = rank(B); % Rank of matrix B
    R_Bs = min(R_B,s);
    K_max = min(q*ceil(R_B/s), n-R_Bs+1);
    K = ceil(n/R_Bs);
    Residue = 1e-18*ones(1, K_max-K+1);
    j=1;
    % The true value of K should be reduced by 1 after returning.
    % Ammendment % (r > del) && % Ammendment
    while (r > del) && (K <= K_max)
        b = xf - (A^K)*x0;
        R = CtrlMatrix(A,B,K);
        u = POMP(R,b,K,s);
        r = norm(xf-(A^K)*x0-R*u);
        Residue(j) = r;
        j = j+1;
        K = K+1;
    end
end