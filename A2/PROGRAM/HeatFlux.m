function [qheatGLO, posgp]= HeatFlux(COOR,CN,TypeElement,ConductMglo,d)
% Function to calculate heat flux in every Gauss point of the domain
% 2D problem is defined so for very element 8 heat fluxes must be obtained.

% Quadrature and shape functions

nnode = size(COOR,1);  % Number of nodes
ndim = size(COOR,2);   % Spatial Dimension of the problem  (2 or 3)
nelem = size(CN,1);   % Number of elements 
nnodeE = size(CN,2) ; %Number of nodes per element 

TypeIntegrand = 'K';

[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,...
    TypeIntegrand) ; 
ngaus = length(weig);
qheatGLO = zeros(ngaus*ndim,nelem); 
auxQ = zeros(2,1);

for e = 1:nelem
    % Define d_el vector
    d_el = zeros(nnodeE, 1);
    XeT = zeros(nnodeE,ndim);
    for j = 1:nnodeE
        d_el(j) = d(CN(e, j), 1);
        XeT(j, :) = COOR(CN(e,j), :);
        Xe = XeT';
    end
    % Make B*d_el products
    for g = 1:ngaus        % for every point of gauss
        BeXi = dershapef(:,:,g) ; 
        % Jacobian Matrix 
        Je = Xe*BeXi' ; 
        % JAcobian 
        % Matrix of derivatives with respect to physical coordinates 
        Be = inv(Je)'*BeXi ; 
        auxQ(1,1) = Be(1,:)*d_el;
        auxQ(2,1) = Be(2,:)*d_el;
        auxQ = -ConductMglo(:,:,e)*auxQ;
        qheatGLO(2*g-1, e) = auxQ(1);
        qheatGLO(2*g, e) = auxQ(2);
    end
end

end

