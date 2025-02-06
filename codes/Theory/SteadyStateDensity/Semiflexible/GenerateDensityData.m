clc; clear;

N = 6; Dr = 5.0; ks = 8.0; kb = 8.0; nu = 1.0; D=1; tau = 1/Dr;

Ms = ChainConnectivityMatrix(N);
Mb = ChainBendingMatrix(N);

alpha = AlphaMiddleActive(N);
file1 = fopen("./../../Data/Semiflexible/Theory/N6-ks8-kb16/alpha_3.dat", 'w'); 

x = linspace(0,100,51);
v = 20*(1+sin(2*pi*x/100));
rho = SteadyStateDensity(N, ks, kb, Ms, Mb, alpha, nu, tau, v, D);


for i = 1:length(x)
    fprintf(file1,"%f \t %f \n", x(i), rho(i));
end
fclose(file1);