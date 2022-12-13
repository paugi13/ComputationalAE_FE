function [X, Def] = ComputeAnalyticalDeformation()
%This funtions gives the Deformation distribution along a beam with a
%square halo cross section subjected to a uniform distributed force.

%TEST DATA
E=70e9;
L=2;
t=-500e3;
hy=0.25;
e=0.05;

N = 100; %Number of data points 
I = 1/12*(hy^4-(hy-2*e)^4); % Beam's Section Inertia

%Vectors decraration
x_beam = zeros(N,1);
def_beam = zeros(N,1);

%Symbolic variables & functions definition
syms x 
M = -t*hy*L*x+t*hy*1/2*x^2+t*hy*L*L/2;

for i=1:numel(x_beam) 
    x_beam(i) = i*L/N;
    expr = M*(x_beam(i)-x);
    res = int(expr,[0 x_beam(i)]);
    def_beam(i) =1/(E*I)*res;
end
X=x_beam; 
Def=def_beam;

%% POSTPROCESSING
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure
hold on
plot(X,Def, 'b');
xlabel('$X$ axis coordinate [$m$]');
ylabel('$Y$ displacements [$m$]');
title('Analytical $Y$ axis displacements');
grid on;
hold off



end