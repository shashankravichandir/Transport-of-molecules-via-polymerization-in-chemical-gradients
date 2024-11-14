clc; clear;

Dr = 5.0; k = 12.0; nu = 1.0; D =1;
gamma_1 = k*nu; tau = 1/Dr;

for Nprime = [2,3,4,5,6,7]
    for str = [1,2,3]
        M = ChainConnectivityMatrix(Nprime);
        switch str
            case 1
                alpha = AlphaMiddleActive(Nprime);
            case 2
                alpha = AlphaAllActive(Nprime);
            case 3
                alpha = AlphaLeadingActive(Nprime);
        end

        [phi, lambda] = eig(M);
        phi = inv(phi);

        for i=1:Nprime
            gamma(i) = gamma_1*lambda(i,i);
        end

        if Nprime==2 && str==1
            N = Nprime;
            k = str;
            eps = EpsilonCalculation(alpha, Nprime, phi, tau, gamma);
        else
            N = [N, Nprime];
            k = [k, str];
            eps = [eps, EpsilonCalculation(alpha, Nprime, phi, tau, gamma)];
        end
    end
end

file1 = fopen("./../Data/Theory/StateDiagram/Dr5k12.dat", 'w'); 
for i = 1:length(N)
    fprintf(file1,"%f \t %f \t %f \n", N(i), k(i), eps(i));
end
fclose(file1);


% names = {'Interior Active', 'All Active', 'Leading Active'};
% fig = figure;
% scatter(N,k, 300, eps, "filled");
% hold on;
% box on;
% set(gca,'ytick',1:3,'yticklabel',names)
% c = colorbar;
% c.Position = [0.85 0.2 0.03 0.7];
% colormap summer;
% xticks(2:6);
% xlim([1.5 7.5]);
% ylim([0.5 3.5]);
% xlabel("N");
% %text(6.6,-0.6,"Antichemotactic", 'FontSize', 12)
% %text(6.7, 2.6,"Chemotactic", 'FontSize', 12)
% 
% set(fig, 'units','points', 'Position', [50, 50, 700, 400]);
% set(gca, 'linewidth', 2, 'OuterPosition', [0 0 1 1], 'Position', [0.25 0.2 0.55 0.7], 'FontSize', 18);
% set(findall(fig,'-property','Interpreter'),'Interpreter','latex');
% set(get(gca, 'YLabel'), 'Rotation', 0, 'VerticalAlignment', 'middle', 'FontSize', 28,'Units', 'Normalized', 'Position', [-0.15,0.5,0]);
% set(get(gca, 'XLabel'), 'Rotation', 0, 'VerticalAlignment', 'middle', 'FontSize', 28,'Units', 'Normalized', 'Position', [0.5,-0.15,0]);
