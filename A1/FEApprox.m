clc;
clear; close all;

n = [5 10 20 40];
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
errorDer = zeros(1, length(n));

for i = 1:length(n)
[errorVector, errorDerVector] = ErrorCalculator(n(i), coord(i, 1:(n(i)+1)), u(i, 1:(n(i)+1)));
error(1, i) = sqrt(max(errorVector));
errorDer(1, i) = sqrt(max(errorDerVector));
end

%% Post-process
% % plot multiple solutions
% figure
% hold on
% plot(coord(1,1:(n(1)+1)), u(1,1:(n(1)+1)));
% plot(coord(2,1:(n(2)+1)), u(2,1:(n(2)+1)));
% plot(coord(3,1:(n(3)+1)), u(3,1:(n(3)+1)));
% plot(coord(4,1:(n(4)+1)), u(4,1:(n(4)+1)));
% xlabel('x (m)');
% ylabel('u (m)');
% title('Finite element approximations');
% legend('nElem = 5', 'nElem = 10', 'nElem = 20', 'nElem = 40', 'location', 'northwest');
% grid on
% hold off

load('exactSolExp.mat');
% comparison with exact solution
figure
subplot(2,2,1);
hold on
fplot(exactSol, [0 1], 'b');
plot(coord(1,1:(n(1)+1)), u(1,1:(n(1)+1)), 'r');
xlabel('x (m)');
ylabel('u (m)');
grid on
title('nEl = 5');
legend('Exact sol.', 'FE sol.', 'location', 'northwest');
hold off

subplot(2,2,2);
hold on
fplot(exactSol, [0 1], 'b');
plot(coord(2,1:(n(2)+1)), u(2,1:(n(2)+1)), 'r');
xlabel('x (m)');
ylabel('u (m)');
grid on
title('nEl = 10');
hold off

subplot(2,2,3);
hold on
fplot(exactSol, [0 1], 'b');
plot(coord(3,1:(n(3)+1)), u(3,1:(n(3)+1)), 'r');
xlabel('x (m)');
ylabel('u (m)');
grid on
title('nEl = 20');
hold off

subplot(2,2,4);
hold on
fplot(exactSol, [0 1], 'b');
plot(coord(4,1:(n(4)+1)), u(4,1:(n(4)+1)), 'r');
xlabel('x (m)');
ylabel('u (m)');
grid on
title('nEl = 40');
hold off






% element size vector

elSizeVector = zeros(length(n), 1);

for i=1:length(elSizeVector)
    elSizeVector(i, 1) = coord(i, 2) - coord(i,1);
end

% plot error vs nElement
figure
plot(log(n), log(error), 'b');
hold on
plot(log(n), log(errorDer), 'r');
% loglog(n, errorDer);
xlabel('log(nElements)');
ylabel('log(\epsilon_{max})');
grid on
title('Total error vs number of elements');
legend('||e||', '||e''||');
hold off


% plot error vs element size
figure
plot(log10(elSizeVector), log10(error), 'b');
hold on
plot(log10(elSizeVector), log10(errorDer), 'r');
xlabel('log(elementSize)');
ylabel('log(\epsilon_{max})');
grid on
title('Total error vs element size');
legend('||e||', '||e''||', 'location', 'northwest');
hold off

