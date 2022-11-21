function Fs = ComputeFs(COOR,CN,TypeElement, fNOD) 
% This subroutine   returns the  heat source contribution (Fs)    to the
% global flux vector. Inputs
% --------------
% 1. Finite element mesh 
% -------------------
% COOR: Coordinate matrix (nnode x ndim)
% CN: Connectivity matrix (nelem x nnodeE)
% TypeElement: Type of finite element (quadrilateral,...)
% -----------
% 2. Vector containing the values of the heat source function at the nodes
% of the mesh
% -----------
%  fNOD (nnode x 1)  %  
%%%%
 
% Dimensions of the problem 
nnode = size(COOR,1); ndim = size(COOR,2); nelem = size(CN,1); nnodeE = size(CN,2) ;     

Fs = zeros(nnode,1) ;  
%......
% warning('You must program the assembly of the flux vector  Fs !!')
     % internal sources. 
TypeIntegrand = 'K';
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand);
XeT = zeros(nnodeE,ndim);

for e=1:nelem  
    for j = 1:nnodeE
        fe = fNOD(CN(e,j));
        XeT(j, :) = COOR(CN(e,j), :);    % transposed matrix
    end
    Xe = XeT';
    Fse = ComputeFseVector(fe,weig,shapef,dershapef,Xe);
    for a = 1:nnodeE
        A = CN(e,a);
        Fs(A) = Fs(A) + Fse(a);
    end
end
