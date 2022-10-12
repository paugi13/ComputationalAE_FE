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

%% Post-process
figure
hold on
plot(coord(1,1:(n(1)+1)), u(1,1:(n(1)+1)));
plot(coord(2,1:(n(2)+1)), u(2,1:(n(2)+1)));
plot(coord(3,1:(n(3)+1)), u(3,1:(n(3)+1)));
plot(coord(4,1:(n(4)+1)), u(4,1:(n(4)+1)));
xlabel('x (m)');
ylabel('u (m)');
hold off



