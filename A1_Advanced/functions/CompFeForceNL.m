function [Fe, eps, sigma] = CompFeForceNL(e, d_k, COOR, CN, AreaFUN, StressFUN)

% Define N for quadratic type. 
elCOOR = [COOR(e); COOR(e+1)];
xi1 = 1/sqrt(3);
xi2 = -1/sqrt(3);
N = zeros(2);
N(1, :) = [1-xi1, 1+xi1];
N(2, :) = [1-xi2, 1+xi2];
N = 0.5*N;

he = elCOOR(2) - elCOOR(1);
B = 1/he*[-1 1];

% meanX = (COOR(e) + COOR(e+1))/2;
eps = B*d_k(CN(e,:));
sigma = StressFUN(eps);

Asum = 0;
for i=1:size(N,1)
    xPh = N(i,:)*elCOOR;
    Asum = Asum + AreaFUN(xPh);
end

% Gauss summation
Asum = he/2*Asum;

% Simple calculations
% meanX = (elCOOR(1)+elCOOR(2))/2;
% A2 = AreaFUN(meanX);
Fe = B.'*Asum*sigma;

% meanX = (elCOOR(1)+elCOOR(2))/2;
% A = AreaFUN(meanX);
% Fe = B.'*A*sigma;


end

