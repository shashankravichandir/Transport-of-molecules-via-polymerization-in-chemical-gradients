function I = integrateMFPT(A, B,a,b)
    integrand = @(x) 2*A(x)./B(x);
    psi = @(x) exp(integral(integrand,a,x,"ArrayValued",1));
    theta = @(x) 1./psi(x);

    integrand = @(x) psi(x)./B(x);
    phi = @(x) integral(integrand,a,x,"ArrayValued",1)./psi(x);

    I = @(x) 2*(integral(theta,a,x,"ArrayValued",1)*integral(phi,x,b,"ArrayValued",1) - integral(theta,x,b,"ArrayValued",1)*integral(phi,a,x,"ArrayValued",1))/integral(theta,a,b,"ArrayValued",1);
end