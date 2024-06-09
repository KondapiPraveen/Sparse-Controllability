% Othogonal Matching Pursuit (OMP)
% Input : A - Measurement Matrix
% Input : y - Measurement Vector
% Input : s - sparsity
% Output : x - sparse vector
function x = OMP(A,y,s)
    m = size(A,2);
    x = zeros(m,1); % Initialization
    V = 1:m; Avl = V;
    Supp = []; % Support Set
    Anorm = vecnorm(A);
    %A = A./Anorm;
    for i=1:s
        [~,ij] = max(abs((A(:,Avl).')*(y-A*x)));
        j = Avl(ij);
        Supp = union(Supp,j);
        Avl = setdiff(Avl,j);
        x(Supp) = pinv(A(:,Supp))*y;
    end
end