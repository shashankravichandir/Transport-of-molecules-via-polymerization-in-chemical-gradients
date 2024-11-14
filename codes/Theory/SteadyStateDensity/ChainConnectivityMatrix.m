function M = ChainConnectivityMatrix(N)
    M = zeros(N,N);
    for i=1:N
        for j=1:N
            if i==j
                if i==1 || i==N
                    M(i,j) = 1;
                else
                    M(i,j) = 2;
                end
            elseif j==i-1 || j==i+1
                M(i,j) = -1;
            end
        end
    end
end