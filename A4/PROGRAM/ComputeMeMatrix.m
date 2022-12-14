function Me = ComputeMeMatrix(densM,weig, shapef, dershapef,Xe)
% Given  % celas: Elasticity Matrix
% weig :   Vector of Gauss weights (1xngaus), % dershapef:   Array with the derivatives of shape functions, with respect to  element coordinates (ndim x nnodeE x ngaus),  Xe: Global coordinates of the nodes of the element,  
% this function returns the element stiffness matrix Ke
ndim = size(Xe,1) ; ngaus = length(weig) ; nnodeE = size(Xe,2)  ;  
Me = zeros(nnodeE*ndim,nnodeE*ndim) ; 
for  g = 1:ngaus
    % Matrix of derivatives for Gauss point "g"
    BeXi = dershapef(:,:,g) ; 
    % Matrix of shape functions at point "g"
    NeSCL = shapef(g,:) ; 
    % Matrix of shape functions at point "g" (for vector-valued fields)
    Ne = StransfN(NeSCL,ndim) ; 
    % Jacobian Matrix 
    Je = Xe*BeXi' ; 
    % JAcobian 
    detJe = det(Je) ;
    NeT = Ne.';
    Me = Me + weig(g)*detJe*(densM*NeT*Ne); 
end