clc;
clear; close all;

g = 0.01;
L = 1;
rho = (pi^2)/(L^2);
s = g*rho^2;
b = g*(pi^2)/L;

% Define the equation
syms x u(x)
f = -s*x^2;

%% FINITE ELEMENT APPROXIMATION
x0 = 0;
x1 = L;
nelem = 10;
nnod = nelem+1;
coord = linspace(x0, x1, nnod);

CN = zeros(nelem, 2);

for i=1:nnod
    CN(i,1) = i;
    CN(i,2) = i+1;
end

w = 1;
xi1 = 1/sqrt(3);
xi2 = -1/sqrt(3);

N1 = [1-xi1, 1+xi1];
N2 = [1-xi2, 1+xi2];

lambda = N1.'*N1*w + N2.'*N2*w;

K = AssemblyK(COOR,CN,lambda, rho);



