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

lambda = N(1, :).'*N(1, :)*w + N(2, :).'*N(2, :)*w;

K = AssemblyK(COOR,CN,lambda, rho);

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

function [Fe] = Compute_Fe_Force(f,N);
    he/2*((1/2*[(1-(1/sqrt(3))) (1+(1/sqrt(3)))])*(-s*x(i)^2)



end



