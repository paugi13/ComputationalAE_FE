clc;
clear; close all;

n = [5 10 15 20 25 30 35 40];
u = zeros(length(n), n(length(n))+1);
coord = zeros(length(n), n(length(n))+1);

for i=1:length(n)
    [u_aux, coord_aux] = SolCalculator(n(i));
    u(i, 1:length(u_aux)) = u_aux;
    coord(i, 1:length(coord_aux)) = coord_aux; 
end

%% Error calculus
load('exactSolExp.mat');
error = zeros(1, length(n));

for i = 1:length(n)
errorVector = ErrorCalculator(n(i), coord(i, 1:(n(i)+1)), u(i, 1:(n(i)+1)));
error(1, i) = sqrt(max(errorVector));
end

%% Post-process
% plot multiple solutions
figure
hold on
plot(coord(1,1:(n(1)+1)), u(1,1:(n(1)+1)));
plot(coord(2,1:(n(2)+1)), u(2,1:(n(2)+1)));
plot(coord(3,1:(n(3)+1)), u(3,1:(n(3)+1)));
plot(coord(4,1:(n(4)+1)), u(4,1:(n(4)+1)));
xlabel('x (m)');
ylabel('u (m)');
title('Finite element approximations');
legend('nElem = 5', 'nElem = 10', 'nElem = 20', 'nElem = 40', 'location', 'northwest');
hold off

% plot error
figure
semilogy(n, error, 'b');
xlabel('Number of elements');
ylabel('log(\epsilon_{max})');
grid on
title('Total error vs number of elements');





