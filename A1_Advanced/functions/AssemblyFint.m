function [Residual,STRAIN,STRESS] = AssemblyFint(COOR,CN,d_k,StressFUN,...
    AreaFUN)

 % d_k zeros(nnode,1) 
%  COOR = linspace(0,L,nnode)' ;
%  CN = [(1:(nnode-1))',(2:nnode)'] ;
nnodeE = 2;

%% K assembly
Ff = zeros(size(COOR, 1), 1);

eps = zeros(size(CN,1), 1);
sigma = zeros(size(CN,1), 1);

for e=1:size(CN, 1)
    [Fe, eps(e, 1), sigma(e, 1)] = CompFeForceNL(e, d_k, COOR, CN, AreaFUN, StressFUN);
    for a = 1:nnodeE
        A = CN(e,a);
        Ff(A) = Ff(A) + Fe(a);
    end
end
Residual = Ff;
STRAIN = eps;
STRESS = sigma;

end

