function [Residual,STRAIN,STRESS] = AssemblyFint(COOR,CN,d_k,StressFUN,...
    AreaFUN)

 % d_k zeros(nnode,1) 
%  COOR = linspace(0,L,nnode)' ;
%  CN = [(1:(nnode-1))',(2:nnode)'] ;
nnodeE = 2;

%% K assembly
Ff = zeros(size(COOR, 2), 1);

for e=1:size(CN, 1)
    [Fe, eps, rho] = CompFeForceNL(e, d_k, COOR, CN, AreaFUN, StressFUN);
    for a = 1:nnodeE
        A = CN(e,a);
        Ff(A) = Ff(A) + Fe(a);
    end
end

Residual = Ff;
STRAIN = eps;
STRESS = rho;

end

