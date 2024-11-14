clc; clear;

N = 6; Dr = 5.0; k = 8.0; nu = 1.0; D=1;
gamma_1 = k*nu; tau = 1/Dr;

M = ChainConnectivityMatrix(N);

alpha = AlphaEndsActive(N);

[phi, lambda] = eig(M);
phi = inv(phi);

for i=1:N
    gamma(i) = gamma_1*lambda(i,i);
end

x = linspace(0,100,51);
v = 20*(1+sin(2*pi*x/100));
rho = SteadyStateDensity(alpha, N, phi, tau, gamma, v, D);

file1 = fopen("./../Data/Theory/alpha_1_6.dat", 'w'); 
for i = 1:length(x)
    fprintf(file1,"%f \t %f \n", x(i), rho(i));
end
fclose(file1);