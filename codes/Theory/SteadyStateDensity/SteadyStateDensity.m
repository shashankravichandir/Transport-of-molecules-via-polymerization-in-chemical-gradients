function rho = SteadyStateDensity(alpha, N, phi, tau, gamma, v, D)
A1 = 0; A2 = 0; A3 = 0;
    %% Calculating A1
for i=1:N
    A1 = A1 + phi(1,i)^2*alpha(i);
end

%% Calculating A2
for i=1:N
    for l = 1:N
        if l~= 1
            A2 = A2+ phi(l,i)^2*alpha(i);
        end
    end
end

%% Calculating A3
for i=1:N
    for j=1:N
        A3 = A3 + (tau*gamma(j)/(1+tau*gamma(j)))*phi(j,i)^2*alpha(i);
    end
end

epsilon = (A1 + A3 - A2)/A1;

rho = (1 + (tau*A1*v.^2)/(2*D)).^(-epsilon/2);
rho = rho./mean(rho);
end