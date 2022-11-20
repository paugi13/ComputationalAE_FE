function K = ComputeK(COOR,CN,TypeElement, ConductMglo)
%%%%
% This subroutine   returns the global conductance matrix K (nnode x nnode)
% Inputs
% --------------
% 1. Finite element mesh
% -------------------
% COOR: Coordinate matrix (nnode x ndim)
% CN: Connectivity matrix (nelem x nnodeE)
% TypeElement: Type of finite element (quadrilateral,...)
% -----------
% 2. Material
% -----------
%  ConductMglo (ndim x ndim x nelem)  % Array of conductivity matrices
%%%%
 if nargin == 0
     load('tmp1.mat')
 end
 
%% COMPLETE THE CODE ....
%warning('You must program the assembly of the conductance matrix K !!')


% Dimensions of the problem
nnode = size(COOR,1);  % Number of nodes
ndim = size(COOR,2);   % Spatial Dimension of the problem  (2 or 3)
nelem = size(CN,1);   % Number of elements 
nnodeE = size(CN,2) ; %Number of nodes per element 

% Determine Gauss weights, shape functions and derivatives  
TypeIntegrand = 'K'; 
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ; 

% Assembly of matrix K
% ----------------
K = sparse(nnode,nnode);
XeT = zeros(nnodeE,ndim);
for e = 1:nelem 
   for j = 1:nnodeE
       XeT(j, :) = COOR(CN(e,j), :);    % transposed matrix
   end
   Xe = XeT';   % in XeT the coordinates from the element's nodes are stored.
   ConductM = ConductMglo(:,:,e);
   Ke = ComputeKeMatrix(ConductM,weig,dershapef,Xe);
   for a = 1:nnodeE
        for b = 1:nnodeE
        % Same working principle as the aerospace structures' one.
            A = CN(e,a);
            B = CN(e,b);
            K(A,B) = K(A,B) + Ke(a,b);
        end
   end
end

