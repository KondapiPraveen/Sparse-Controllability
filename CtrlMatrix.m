function R = CtrlMatrix(A,B,t)
    [n,m] = size(B);
    R = zeros(n,m*t);
    for i=0:t-1
        R(:,1+i*m:(i+1)*m) = A^(i)*B;
    end
end