function [qheatGLO, posgp]= HeatFlux(COOR,CN,TypeElement,ConductMglo,d)
% Function to calculate heat flux in every Gauss point of the domain
% 2D problem is defined so for very element 8 heat fluxes must be obtained.

% Quadrature and shape functions
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,...
    TypeIntegrand) ; 

nnode = size(COOR,1);  % Number of nodes
ndim = size(COOR,2);   % Spatial Dimension of the problem  (2 or 3)
nelem = size(CN,1);   % Number of elements 
nnodeE = size(CN,2) ; %Number of nodes per element 

qheatGLO = zeros(ngaus*ndim,nelem); 

for i = 1:nelem
    % Define d_el vector
    d_el = zeros(nnodE, 1);
    for j = 1:nnodeE
        d_el(j) = d(CN(i, j), 1);
    end
    % Make B*d_el products
    for g = 1:nnodeE
        qheatGLO(2*g-1, i) = dershapef(1,:,g)*d_el;
        qheatGLO(2*g, i) = dershapef(2,:,g)*d_el;
    end
end

end

