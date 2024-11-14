function tm = MeanFirstPassageTime(N, Dr, k, nu, D, alpha)
    a = 25; b = 125;
    gamma_1 = k*nu; tau = 1/Dr;
    M = ChainConnectivityMatrix(N);

    [phi, lambda] = eig(M);
    phi = inv(phi);
    for i=1:N
        gamma(i) = gamma_1*lambda(i,i);
    end

    [eps,S2] = epsilon(alpha, N, phi, tau, gamma);

    v = @(x) 20*(1+sin((2*pi*x)/100));
    dv = @(x) (2*pi/5)*cos((2*pi*x)/100);

    %% Calculating psi
    Diff = @(x) (1/N)*(D + (tau/2)*S2*v(x).^2);
    dDiff = @(x) (1/N)*((tau/2)*S2*2*v(x).*dv(x));

    A = @(x) ((1-eps/2)*dDiff(x));
    B = @(x) (2*Diff(x));
    T = integrateMFPT(A,B,a,b);
    tm = T(50);
end