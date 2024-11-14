function [eps,S2] = epsilon(alpha, N, phi, tau, gamma)
S1 = 0; S2 = 0; S3 = 0;
    %% Calculating S2
for i=1:N
    S2 = S2 + phi(1,i)^2*alpha(i);
end

%% Calculating S1
for i=1:N
    for l = 1:N
        if l~= 1
            S1 = S1+ phi(l,i)^2*alpha(i);
        end
    end
end

%% Calculating S3
for i=1:N
    for j=1:N
        S3 = S3 + (tau*gamma(j)/(1+tau*gamma(j)))*phi(j,i)^2*alpha(i);
    end
end

eps = (S2 + S3 - S1)/S2;
end