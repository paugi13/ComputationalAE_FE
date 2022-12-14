% Simple program to plot displacements along one edge
clc;
clear;
close all;

m1 = [0 0
1 -0.00086799997
2 -0.0024580001];

m2 = [0 0
0.2 -0.000239
0.40000001 -0.00080500002
0.60000002 -0.001642
0.80000001 -0.0026809999
1 -0.0038689999
1.2 -0.0051589999
1.4 -0.006513
1.6 -0.0079009999
1.8 -0.0093019996
2 -0.010703];

m3 = [0 0
0.133333 -0.000137
0.26666701 -0.00043700001
0.40000001 -0.00089600001
0.533333 -0.001482
0.66666698 -0.002177
0.80000001 -0.002961
0.93333298 -0.003817
1.0666699 -0.0047280001
1.2 -0.005682
1.33333 -0.0066669998
1.46667 -0.0076720002
1.6 -0.008688
1.73333 -0.0097110001
1.86667 -0.010734
2 -0.011756];

m4 = [0 0
0.1 -9.4000003e-05
0.2 -0.00028099999
0.30000001 -0.00056999997
0.40000001 -0.00094300002
0.5 -0.001392
0.60000002 -0.001906
0.69999999 -0.002478
0.80000001 -0.0030990001
0.89999998 -0.003762
1 -0.0044590002
1.1 -0.0051839999
1.2 -0.0059329998
1.3 -0.006699
1.4 -0.0074769999
1.5 -0.0082649998
1.6 -0.0090589998
1.7 -0.0098559996
1.8 -0.010654
1.9 -0.011451
2 -0.012247];

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

fig1 = figure(1);
hold on
plot(m1(:,1), m1(:,2));
plot(m2(:,1), m2(:,2));
plot(m3(:,1), m3(:,2));
plot(m4(:,1), m4(:,2));
grid on
xlabel('$X$ axis coordinate [$m$]');
ylabel('$Y$ displacements [$m$]');
title('$Y$ axis displacements along upper right edge');
legend('$n_{ex} = 2$ $d = 0,1~m$', '$n_{ex} = 10$ $d = 0,05~m$', '$n_{ex} = 15$ $d = 0,05~m$', ...
    '$n_{ex} = 20$ $d = 0,025~m$', 'location', 'southwest');
hold off