% Selecting LI Columns by minimizing Tr(W_S^{-1})
% Input -  
% System Matrices A, B(Rank(B)=n)
% Control Time : K >= n/s
% Sparsity Level : s
% Output - S (Binary Matrix)
function [S,S_k] = FullBLI(A,B,K,s)
    % Initialization
    m = size(B,2); n = size(A,1);
    S_k = zeros(K,m); e_t = 1e-15;
    S = []; IW_S = (1/e_t)*eye(n);
    rankC = zeros(K,1); % Collection of rank information A^iB
    R = zeros(n,m*K);
    % Find rank(A^iB)
    I = B;
    for i=0:K-1
        R(:,i*m+1:(i+1)*m)=I;
        rankC(i+1) = rank(I);
        if i<K-1
            I = A*B;
        end
    end
    % Adding LI Columns that minimize Average Energy
    NumLI = 0; % Count of # LI columns in each iteration
    for i = K-1:-1:0
        V = 1:m; % The set of available actuators
        AiB = R(:,i*m+1:(i+1)*m); % The matrix A^iB \in (n,m)
        l_i = min(s,rankC(i+1)-NumLI); % Max # LI columns that can be added
        for z = 1:l_i
            v = AiB(:,V);
            X = v.'*IW_S;
            Y = sum(X.*(v.'),2);
            UpVec = (vecnorm(X,2,2).^2)./(1+Y); % Update Vector
            [~, l] = max(UpVec);
            p = V(l); % actuator no. at that time instant
            V = setdiff(V,p); % Remove the selected actuator
            % k = ceil(p/m); j = p-(k-1)*m; % Find the time step and actuator no.
            j = i*m+p; % Actual actuator no. in grand scale
            if sum(S_k(i+1,:)) < s % Check if the new actutator can be added
                S_k(i+1,p)=1;
                S = union(S,j);
                v1 = AiB(:,p);
                IW_S = IW_S - ((IW_S*(v1*v1.')*IW_S)/(1+v1.'*IW_S*v1));
            end
            %{
            fprintf(strcat(num2str(p),' \n'));
            str = sprintf("Current Rank r = %d\n",rank(R(:,S)));
            fprintf(str);
            %}
        end
        NumLI = NumLI+l_i;
        %fprintf('--------------------------------\n');
    end
    if rank(R(:,S))<n
        fprintf("Rank Deficient \n");
    end
    %{
    str = sprintf("The trace value is %.3f\n",trace(IW_S));
    fprintf(str);
    %}
end