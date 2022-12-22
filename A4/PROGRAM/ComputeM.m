function M = ComputeM(COOR,CN,TypeElement,densglo) 
% THis subroutine assembles the global mass matrix of the structre
if nargin == 0
    load('tmp1.mat');
end
nnode = size(COOR,1); 
ndim = size(COOR,2); 
nelem = size(CN,1); 
nnodeE = size(CN,2) ;  
% nstrain = size(celasglo,1) ;
% Shape function routines (for calculating shape functions and derivatives)
TypeIntegrand = 'K';
[weig,~,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;
% Assembly of matrix K
% ----------------
M = sparse([],[],[],nnode*ndim,nnode*ndim,nnodeE*ndim*nelem) ;
for e = 1:nelem
    dens = densglo(e) ;  % Density of element "e"
    CNloc = CN(e,:) ;   % Coordinates of the nodes of element "e"
    Xe = COOR(CNloc,:)' ;     % Computation of elemental mass matrix
    Me = ComputeMeMatrix(dens,weig, shapef, dershapef,Xe);
    
    for anod = 1:nnodeE
        a = Nod2DOF(anod,ndim) ; 
            Anod = CN(e,anod) ; 
            A = Nod2DOF(Anod,ndim) ;         
        for bnod = 1:nnodeE
            b = Nod2DOF(bnod,ndim) ;
            Bnod = CN(e,bnod);
            B = Nod2DOF(Bnod,ndim) ;
            M(A,B) = M(A,B) + Me(a,b);
        end
    end
       % -----------------------------------------------------------------------
    if mod(e,10)==0  % To display on the screen the number of element being assembled
        disp(['e=',num2str(e)])
    end
end


