clc;
clear;
close all;

% Add subroutines to path
addpath('ROUTINES_AUX');

% Load analysis data
name = 'P4ModesCase2'; 
load(name);

% Calculate first 25 modes
neig = 25;

npt3 = size(K, 1);  % number of dofs
% Changes in noation
M = M(DOFl, DOFl);
K = K(DOFl, DOFl);
d = d(DOFl);
[MODES, FREQ] = UndampedFREQ(M, K, neig);
FREQhz = FREQ/(2*pi);

% GID postprocess
GidPostProcessModes(COOR,CN,TypeElement,MODES,posgp,NameFileMesh,DOFl);

%% Dynamic response analysis
% Initial velocity condition
dDer = zeros(length(DOFl), 1);
xi = 0;      % Damping factor
m = 40;
switch name
    case 'P4ModesCase1'
        T1 = 2*pi/FREQ(1);
        totalT = m*T1;
    case 'P4ModesCase2'
        T1 = 2*pi/FREQ(5);
        totalT = m*T1;
end

nstep = 500;
t = linspace(0, totalT, nstep);
t = t';

dTimeVector = zeros(length(DOFl), nstep);
qi0Vector = zeros(1, neig);
for j = 1:size(t, 1)
    disp(['time step = ', num2str(j)]);
    for i = 1:neig
        qi0 = MODES(:, i)'*M*d;
        if j == 1
            qi0Vector(i) = qi0;
        end
        qi0Der = MODES(:, i)'*M*dDer;
        freqBar = FREQ(i)*sqrt(1-xi^2);
        % First parenthesis
        f1 = qi0*cos(freqBar*t(j)) + ((qi0Der + xi*FREQ(i)*qi0)/freqBar)*sin(freqBar*t(j));
        fT = exp(-xi*FREQ(i)*t(j))*f1;
        dTimeVector(:, j) = dTimeVector(:, j) + MODES(:, i)*fT;
    end
end

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Visualize mode amplitudes
qi0Vector = abs(qi0Vector);
fig1 = figure(1);
hold on
bar(qi0Vector);
xlabel('Mode');
ylabel('Amplitude $q_i^0$');
grid on
hold off

% dTimeVector -> displacements of free nodes.
% total matrix including restricted nodes must be assembled.
DISP = zeros(npt3, nstep);
DISP(DOFl, :) = dTimeVector;
NAME_INPUT_DATA = 'DYN-sol2';
GidPostProcessDynamic(COOR,CN,TypeElement,DISP,NAME_INPUT_DATA,posgp,NameFileMesh,t);
