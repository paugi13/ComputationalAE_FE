clear; clc;
close all;

%% EXACT SOLUTION OF THE DIFFERENTIAL EQUATION.

g = 0.01;
L = 1;
rho = (pi^2)/(L^2);
s = g*rho^2;
b = g*(pi^2)/L;

% Define the equation
syms x u(x)
f = -s*x^2;
eq = diff(u,x,2) + rho*u == -f;

% Boundary conditions
Du = diff(u,x);
c1 = u(0) == -g;
c2 = Du(L) == b;

%Solve the equation
exactSol = dsolve(eq, c1, c2);
exactSol = vpa(exactSol,4);
% exactSol = 0.01*cos(pi*x) + 0.098696*x^2 + pi/100*sin(pi*x) - 0.02
figure
fplot(exactSol, [0 1], 'b');
xlabel('x (m)');
ylabel('u (m)');
title('Axial displacements along \Omega');
grid on

%% Derivation in terms of generic functions
N = [1,x] ; % Basis functions (polynomials)
u1 = GalerkinMethod(N,f,b,g,L,rho);
N = [1,x,x^2];
u2 = GalerkinMethod(N,f,b,g,L,rho);
N = [1,x,x^2,x^3];
u3 = GalerkinMethod(N,f,b,g,L,rho);
N = [1,x,x^2,x^3,x^4];
u4 = GalerkinMethod(N,f,b,g,L,rho);

figure
fplot(exactSol, [0 1], 'b');
hold on
xlabel('x (m)');
ylabel('u (m)');
fplot(u1, [0 1]);
fplot(u2, [0 1]);
fplot(u3, [0 1]);
fplot(u4, [0 1]);
%, 'color', [0.3 0.9 0.5]
legend('Exact solution', 'N = [1,x]', 'N = [1,x,x^2]', 'N = [1,x,x^2,x^3]'...
    , 'N = [1,x,x^2,x^3,x^4]', 'location', 'northwest');
title('Displacements using matrix forms of the Principle of Virtual Work');

hold off


