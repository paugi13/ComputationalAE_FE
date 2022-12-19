clc;
clear;
close all;

% Add subroutines to path
addpath('ROUTINES_AUX');

% Load analysis data
load('P4ModesCase1');

neig = 25;

Mll = M(DOFl, DOFl);
Kll = K(DOFl, DOFl);
[MODES, FREQ] = UndampedFREQ(Mll, Kll, neig);

GidPostProcessModes(COOR,CN,TypeElement,MODES,posgp,NameFileMesh,DOFl);