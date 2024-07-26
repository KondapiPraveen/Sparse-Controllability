function MA = RandomCluster(p,NSys)
    n = 2*p;
    MA = zeros(n,n,NSys);
    Aii = zeros(2,2,NSys);
    range1=0.25; bias = 0;
    range2=0.25; V = 1:n;
    range3 = 0.005;
    for i=1:p
        Aii(1,1,:) = rand(1,1,NSys)*range1 + bias;
        Aii(2,2,:) = rand(1,1,NSys)*range1 + bias;
        Aii(1,2,:) = rand(1,1,NSys)*range2;
        Aii(2,1,:) = rand(1,1,NSys)*range2;
        BlckEtr = 2*(i-1)+1:2*i; OffDiag = [1:2*(i-1), 2*i+1:n];
        MA(BlckEtr,BlckEtr,:) = Aii;
        MA(BlckEtr,OffDiag,:) = rand(2,2*(p-1),NSys)*range3;
    end
end