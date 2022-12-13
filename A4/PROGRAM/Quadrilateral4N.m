function [Ne, BeXi] = Quadrilateral4N(xi,eta)
% Quadrilateral element matrix of shape functions. 
Ne = 1/4*[(1-xi)*(1-eta) (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta)];
BeXi = 1/4*[-(1-eta) (1-eta) (1+eta) -(1+eta)
            -(1-xi) -(1+xi) (1+xi) (1-xi)];
end

