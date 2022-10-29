function [Fe, eps, rho] = CompFeForceNL(e, d_k, COOR, CN, AreaFUN, StressFUN)

B = 0.5*[-1 1];
meanX = (COOR(e+1) + COOR(e));
eps = B*d_k(CN(e,:));

A = AreaFUN(meanX);
rho = StressFUN(eps);
Fe = B.'*A*rho*2;    % Fe 2x1


end

