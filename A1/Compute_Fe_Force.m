function Fe = Compute_Fe_Force(f,N,e,coord)
% Function to compute element forces Fe

w = 1;
wcoord = [coord(e); coord(e+1)];
he = wcoord(2) - wcoord(1);
N = 0.5*N;

f = zeros(size(N,1),1);

for i=1:size(N,1)
    xPh = N(i,:)*wcoord;
    f(i,1) = subst(f, x, xPh);
end

% xPh1 = N(1,:)*wcoord;
% f1 = subst(f, x, xPh1);
% xPh2 = N(1,:)*wcoord;
% f2 = subst(f, x, xPh2);

Fe = he/2*(w*N(1,:).'*f(1)+w*N(2,:).'*f(2));
end

