function [u] = FE_Resolution(Ff,K,g,nnod)
r = 1;
l = 2:length(nnod);
dr = -g;
dl = K(l,l)\(Ff(l)+K(l,r)*g);
u = [dr,dl];
end