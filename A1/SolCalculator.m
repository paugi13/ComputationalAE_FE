function [u, coord] = SolCalculator(nelem)
% Function that returns 1D solution given a number of elements 


%% Numerical data

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
nnod = nelem+1;
coord = linspace(x0, x1, nnod);
he = coord(2)-coord(1);
CN = zeros(nelem, 2);

% Connectivity matrix
for i=1:nelem
    CN(i,1) = i;
    CN(i,2) = i+1;
end

w = 1;
xi1 = 1/sqrt(3);
xi2 = -1/sqrt(3);

%% K assembly
N = zeros(2);
N(1, :) = [1-xi1, 1+xi1];
N(2, :) = [1-xi2, 1+xi2];
N = 0.5*N;

lambda = he/2*(N(1, :).'*N(1, :)*w + N(2, :).'*N(2, :)*w);

K = AssemblyK(coord,CN,lambda, rho);

%% F assembly

nnodeE = 2;
Ff = zeros(nnod,1);

for e=1:nelem
    Fe = Compute_Fe_Force(f,N,e,coord);
    for a = 1:nnodeE
        A = CN(e,a);
        Ff(A) = Ff(A) + Fe(a);
    end
end


% Max(xi) = 1
N_1 = [0 1];
Ff([size(Ff,1)-1 size(Ff,1)]) = Ff([size(Ff,1)-1 size(Ff,1)]) + N_1.'*b;

F = Ff;

% Final solution
u = FE_Resolution(F,K,g,nnod);

end

