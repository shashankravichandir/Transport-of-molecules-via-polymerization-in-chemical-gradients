function alpha = AlphaEndsActive(N)
    alpha = zeros(1,N);
    alpha(1) = 1; alpha(N) = 1;
end