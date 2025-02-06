function M = ChainBendingMatrix(N)
    M = zeros(N,N);
    for i=1:N
        for j=1:N
            if i==j
                if i==1 || i==N
                    M(i,j) = 1;
                elseif i==2 || i==N-1
                    M(i,j) = 5;
                else
                    M(i,j) = 6;
                end
            elseif j==i-1
                if i==2||i==N
                    M(i,j) = -2;
                else
                    M(i,j) = -4;
                end
            elseif j==i+1
                if i==1||i==N-1
                    M(i,j) = -2;
                else
                    M(i,j) = -4;
                end
            elseif j==i-2 || j==i+2
                M(i,j) = 1;
            end
        end
    end
end