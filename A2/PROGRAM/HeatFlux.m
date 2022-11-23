function [qheatGLO, posgp]= HeatFlux(COOR,CN,TypeElement,ConductMglo,d)
% Function to calculate heat flux in every Gauss point of the domain
% 2D problem is defined so for very element 8 heat fluxes must be obtained.

% Quadrature and shape functions

ndim = size(COOR,2);   % Spatial Dimension of the problem  (2 or 3)
nelem = size(CN,1);   % Number of elements 
nnodeE = size(CN,2) ; %Number of nodes per element 

TypeIntegrand = 'K';

[weig,posgp,~,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,...
    TypeIntegrand) ; 
ngaus = length(weig);
qheatGLO = zeros(ngaus*ndim,nelem); 
% auxQ = zeros(2,1);

for e = 1:nelem
    % Define d_el vector
    CNloc = CN(e,:);
    Xe = COOR(CNloc, :)';
    d_el = d(CNloc);
    for g = 1:ngaus        % for every point of gauss
        BeXi = dershapef(:,:,g) ; 
        % Jacobian Matrix 
        Je = Xe*BeXi' ; 
        % Matrix of derivatives with respect to physical coordinates 
        Be = inv(Je)'*BeXi ; 
        auxQ = Be*d_el;
        auxQ = -ConductMglo(:,:,e)*auxQ;
        qheatGLO([2*g-1 2*g], e) = auxQ;
    end
end

end
