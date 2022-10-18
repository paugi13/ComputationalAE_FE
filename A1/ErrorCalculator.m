function [errorVector, errorDerVector] = ErrorCalculator(nelem, coord, d)

%% Element function
w = 1;
xi1 = 1/sqrt(3);
xi2 = -1/sqrt(3);

N = zeros(2);
N(1, :) = [1-xi1, 1+xi1];
N(2, :) = [1-xi2, 1+xi2];
N = 0.5*N;

errorVector = zeros(nelem,1);
errorDerVector = zeros(nelem,1);
fPh = zeros(size(N,1),1);

% Load exact symbolic solution
syms x
load('exactSolExp.mat');
B = 0.5*[-1, 1];
% error calculator
for i=1:nelem
    wcoord = [coord(1,i); coord(1,(i+1))];
    he = wcoord(2) - wcoord(1);
    sum = 0;
    de = [d(1,i); d(1,(i+1))];
    for j=1:size(N,1)
        xPh = N(j,:)*wcoord;
        fPh(j,1) = subs(exactSol, x, xPh);
        sum = sum + (w*(fPh(j,1)-B*de))^2;
    end
    errorVector(i,1) = he/2*sum;
end
   
% error's derivative

derExactSol = diff(exactSol, x);

for i=1:nelem
    wcoord = [coord(1,i); coord(1,(i+1))];
    he = wcoord(2) - wcoord(1);
    sum = 0;
    de = [d(1,i); d(1,(i+1))];
    for j=1:size(N,1)
        xPh = N(j,:)*wcoord;
        fPh(j,1) = subs(derExactSol, x, xPh);
        sum = sum + (w*(fPh(j,1)-N(j,:)*de))^2;
    end
    errorDerVector(i,1) = he/2*sum;
end


end

