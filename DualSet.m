% Chooses actuators with Replacement (Weighted) (Adapted to Piecewise Sparsity)
function [c, Supp] = DualSet(V,U,s,ts)
    % V, U are Set of Vectors such that VV.' = I = UU.'
    % s : sparsity level
    % ts : # Time Steps
    % k : # total no. of active actuators
    k = s*ts;
    [n, t] = size(V);
    [l, ~] = size(U);
    m = round(t/ts);
    c = zeros(t,1); A_l = zeros(n); A_u = zeros(l);
    d_l = 1; d_u = (1+sqrt(l/k))/(1-sqrt(n/k));
    Val = 1:t; % Valid Support
    S = zeros(ts,m);
    Supp = []; % Schedule Support
    i=0;
    cst1 = sqrt(k*n);
    cst2 = sqrt(k*l);
    while i < k
        mu_l = i-cst1;
        mu_u = d_u*(i+cst2);
        % Find j
        E_l = real(eig(A_l)); E_u = real(eig(A_u));
        A_li = inv(A_l - (mu_l+d_l)*eye(n));
        A_ui = inv((mu_u+d_u)*eye(l)-A_u);
        A_li2 = A_li^2;
        A_ui2 = A_ui^2;
        if max(imag(A_li(:)),'all') < eps
            A_li = real(A_li);
        end
        
        if max(imag(A_li(:)),'all') < eps
            A_ui = real(A_ui);
        end
        V1 = V(:,Val);
        U1 = U(:,Val);
        cli = CL(V1,d_l,A_li,A_li2,E_l,mu_l);
        cui = CU(U1,d_u,A_ui,A_ui2,E_u,mu_u);
        dif = cli-cui;
        [mdif,q] = max(dif);
        if mdif < 0
            fprintf('Error at DualSet \n');
            break;
        end
        cl = cli(q); cu = cui(q);
        j = Val(q);
        
        %{
        for idx=Val
            cli = CL(V(:,idx),d_l,A_li,A_li2,E_l,mu_l);
            cui = CU(U(:,idx),d_u,A_ui,A_ui2,E_u,mu_u);
            dif = cli-cui;
            if dif > Opt
                j = idx;
                Opt = dif;
                cl = cli; cu = cui;
            end
        end
        %}
        rw = ceil(j/m); p = j-(rw-1)*m;
        % Val = setdiff(Val,j);
        
        if ((sum(S(rw,:)) < s) || (sum(S(rw,:))<=s && S(rw,p) == 1))
            if S(rw,p) == 0
                S(rw,p) = 1;
                Supp = union(Supp,j);
            end
            del = 2/(cl+cu);
            %{
            if del < 0
                fprintf("Sparsity : %d , del : %.2f \n",s,del)
            end
            %}
            ej = zeros(t,1); ej(j) = 1;
            c = c+del*ej;
            i=i+1;
            % Update A
            A_l = A_l + del*V(:,j)*V(:,j).';
            A_u = A_u + del*U(:,j)*U(:,j).';
        end
        if (sum(S(rw,:)) == s)
            fbd = (rw-1)*m+(1:m);
            Val = setdiff(Val,setdiff(fbd,Supp));
        end
    end
    c = (1-sqrt(n/k))*c/k;
end

% Barrier or Potentional Functions
function dv = Phi_L(l, E)
    % E = eig(A);
    dv = sum(1./(E-l));
end

function dv = Phi_U(u, E)
    % E = eig(A);
    dv = sum(1./(u-E));
end

% Parameter Functions
function L = CL(v,d,A,A2,E,mu)
    % L = (v.'*(A2)*v)/(Phi_L(mu+d,E)-Phi_L(mu,E)) - v.'*(A)*v;
    dnm = (Phi_L(mu+d,E)-Phi_L(mu,E));
    X = sum((v.'*A2).*(v.'),2);
    Y = sum((v.'*A).*(v.'),2);
    L = X/dnm - Y;
end

function U = CU(u,d,A,A2,E,mu)
    % U = u.'*(A)*u + (u.'*(A2)*u)/(Phi_U(mu,E)-Phi_U(mu+d,E));
    dnm = (Phi_U(mu,E)-Phi_U(mu+d,E));
    X = sum((u.'*A2).*(u.'),2);
    Y = sum((u.'*A).*(u.'),2);
    U = X/dnm + Y;
end
