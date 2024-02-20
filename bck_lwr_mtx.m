% Toeplitz Block Lower Traingular matrix
% Input : System Matrices - A, B
% Input : Regularizer - R (u.'*R*R*u)
% Input : # Blocks - N
%{
Output : M2
B           0   ...   0
AB          0   ...   0
.
.
A^{N-1}B A^{N-2}B ... B

            R

Output : M1
A A^2 ... A^N
%}
function [M1, M2] = bck_lwr_mtx(A,B,R,N)
    [n,m] = size(B); % State, Input Vector Dimensions
    
    M1 = zeros(n*N,n);
    Z = zeros(n*N,m); % [B AB ... A^{N-1}B]
    for i=0:N-1
        M1(i*n+1:(i+1)*n,:) = A^(i+1);
        if i==0
            Ai = eye(n);
        else
            Ai = M1((i-1)*n+1:i*n,:);
        end
        Z(i*n+1:(i+1)*n,:) = Ai*B;
    end
    
    M2 = zeros((n+m)*N,m*N);
    for i=0:N-1
        M2(i*n+1:n*N,i*m+1:(i+1)*m) = Z(1:end-i*n,:);
    end
    M2(n*N+1:end,:) = R;
end