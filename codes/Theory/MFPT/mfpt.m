clc; clear;

location = "./../Data/Theory/";

Dr = 5.0; k = 8.0; nu = 1.0; D = 1;

for N = [2,3,4,5,6,8,10,12,14,16,18,25,50,100]
    alpha = AlphaAllActive(N);
    if N==2
        tm = MeanFirstPassageTime(N, Dr, k, nu, D, alpha);
    else
        tm = [tm, MeanFirstPassageTime(N, Dr, k, nu, D, alpha)];
    end
end
%%
N  = [2,3,4,5,6,8,10,12,14,16,18,25,50,100];
%%

% N = [2,3,4,5,6,8,10,12,14,16,18];
% file = fopen(location + "mid_active.dat", 'wt');
% l = length(N);
% for i=1:l
%        fprintf(file, '%f\t%f\n', N(i), tm(i));
% end
%file.fclose();





