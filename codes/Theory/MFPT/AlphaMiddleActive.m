function alpha = AlphaMiddleActive(N)
    alpha = zeros(1,N);
    alpha(round(N/2)) = 1;
end