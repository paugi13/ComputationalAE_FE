function [weig,posgp,shapef,dershapef] = Hexahedra8NInPoints

weig = [1 1 1 1 1 1 1 1];
posgp = 1/sqrt(3)*[-1 1 1 -1 -1 1 1 -1
                -1 -1 1 1 -1 -1 1 1
                -1 -1 -1 -1 1 1 1 1];
ndim = 3;
nnodeE = 8;
ngaus = length(weig);
shapef = zeros(ngaus,nnodeE);
dershapef = zeros(ndim,nnodeE,ngaus);
    for g = 1:length(weig)
        xi = posgp(1,g); 
        eta = posgp(2,g);
        zeta = posgp(3,g);
        [Ne, BeXi] = Hexaedra8N(xi,eta,zeta);   % N and B matrices. 
        shapef(g,:) = Ne;           % [Ne functions along columns]
        dershapef(:,:,g) = BeXi;    % BeXi matrix along third dimensions
    end
end

