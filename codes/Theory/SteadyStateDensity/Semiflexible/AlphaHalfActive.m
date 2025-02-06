function alpha = AlphaHalfActive(N)
    alpha = zeros(1,N);
    for i = 1:N/2
        alpha(i) = 1;
    end
end