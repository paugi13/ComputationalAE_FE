function [weig,posgp,shapef,dershapef] = Quadrilateral4NInPoints 
% Four integration points
weig = [1 1 1 1];   % length(weig) can define number of gauss points.
posgp = 1/sqrt(3)*[-1 +1 +1 -1
                    -1 -1 +1 +1];   % first row: xi. second row: eta. 
ndim = 2; nnodeE = 4;
ngaus = length(weig);
shapef = zeros(ngaus,nnodeE);
dershapef = zeros(ndim,nnodeE,ngaus);
    for g = 1:length(weig)
        xi = posgp(1,g); 
        eta = posgp(2,g);
        [Ne, BeXi] = Quadrilateral4N(xi,eta);   % N and B matrices. 
        shapef(g,:) = Ne;           % [Ne functions along columns]
        dershapef(:,:,g) = BeXi;    % BeXi matrix along third dimensions
    end
end

% shapedef contains Ne vectors. 
% dershapedef contains Be matrices. 