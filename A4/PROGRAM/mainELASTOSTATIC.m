clc;
clear; 
% Finite Element Program for Elastostatic problems  
% ECA.
% Technical University of Catalonia
% JoaquIn A. Hdez, October 23-th, 2015
% ---------------------------------------------------
if exist('ElemBnd')==0
    addpath('ROUTINES_AUX') ,
end

 
%%% INPUT  %%% 
% Input data file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NAME_INPUT_DATA = 'BEAM3D' ;  % Name of the mesh file 
%------------------------------------------------------

% PREPROCESS  
[COOR,CN,TypeElement,TypeElementB, celasglo, DOFr,dR,...  
    Tnod,CNb,fNOD,Fpnt,NameFileMesh,typePROBLEM,celasgloINV,DATA, densglo] = ReadInputDataFile(NAME_INPUT_DATA)  ; 

% SOLVER 
% --------------------------------------------
[d, strainGLO, stressGLO,  React, posgp, K, M, DOFl]= SolveElastFE(COOR,CN,TypeElement,TypeElementB, celasglo,  DOFr,dR,...  
    Tnod,CNb,fNOD,Fpnt,typePROBLEM,celasgloINV,DATA, densglo)  ; 

save('P4ModesCase1', 'COOR', 'CN', 'posgp','TypeElement', 'K', 'M', 'NameFileMesh',...
    'd', 'DOFl');

% POSTPROCESS
% --------------------------------------------
% GidPostProcessModes(COOR,CN,TypeElement,MODES,posgp,NameFileMesh,DOFl) 