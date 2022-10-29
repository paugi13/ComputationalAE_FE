function [Fe, eps, sigma] = CompFeForceNL(e, d_k, COOR, CN, AreaFUN, StressFUN)

% Define N for quadratic type. 
xi1 = 1/sqrt(3);
xi2 = -1/sqrt(3);
N = zeros(2);
N(1, :) = [1-xi1, 1+xi1];
N(2, :) = [1-xi2, 1+xi2];
N = 0.5*N;
B = 0.5*[-1 1];

elCOOR = [COOR(e) COOR(e+1)];
% meanX = (COOR(e) + COOR(e+1))/2;
eps = B*d_k(CN(e,:));
sigma = StressFUN(eps);

sum = 0;
for i=1:size(N,1)
    xPh = N(i,:)*elCOOR;
    A = AreaFUN(xPh);
    sum = sum + B.'*A*sigma;    % 2x1
end
Fe = sum;

end

