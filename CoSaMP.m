% Input Variables
% b : Measurement Vector; A : Sensing Matrix
% s : Sparsity level
% max_itr : maximum iteration of the CoSaMP
% tol : termination tolerance of consequetive residues
% Output Variable
% x : input vector with sparsity 's'
function x = CoSaMP(b,A,s,max_itr,tol)
    A_norm = A./vecnorm(A); % col. normalized A
    b_norm = b/norm(b);
    n = size(A,2); % size of the (input) vector 
    x = zeros(n,1);
    itr = 1;
    Supp = []; res_new = b; itr_cond = 1;
    while itr_cond && itr <= max_itr
        res = res_new; % residue
        [~,Ind_2s] = sort(abs(A_norm'*res),'descend'); % largest 2s abs entries
        iSupp = union(Supp, Ind_2s(1:2*s)); % Intermediate Support with 3s
        x_new = zeros(n,1);
        A_iSupp = A(:,iSupp);
        x_new(iSupp) = pinv(A_iSupp)*b;
        [~,Ind_3s] = sort(abs(x_new));
        x_new(Ind_3s(s+1:end)) = 0;
        res_new = b-A*x_new;
        itr_cond = abs(res_new-res)/b_norm > tol;
        Supp = Ind_3s(1:s);
        x = x_new;
        itr = itr+1;
    end
end