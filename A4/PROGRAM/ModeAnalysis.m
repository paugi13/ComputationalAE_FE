clc;
clear;
close all;

% Add subroutines to path
addpath('ROUTINES_AUX');

% Load analysis data
load('P4ModesCase1');

neig = 25;

[MODES, FREQ] = UndampedFREQ(M, K, neig);

GidPostProcessModes(COOR,CN,TypeElement,MODES,posgp,NameFileMesh,DOFl);