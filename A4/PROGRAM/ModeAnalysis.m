clc;
clear;
close all;

% Add subroutines to path
addpath('ROUTINES_AUX');

% Load analysis data
load('P4ModesCase1');

MODES = 3;

GidPostProcessModes(COOR,CN,TypeElement,MODES,posgp,NameFileMesh,DOFl)