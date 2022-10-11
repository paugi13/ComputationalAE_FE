function Fe = Compute_Fe_Force(f,N,e,coord)
% Function to compute element forces Fe

w = 1;
wcoord = [coord(e); coord(e+1)];
he = wcoord(2) - wcoord(1);
N = 0.5*N;

fPh = zeros(size(N,1),1);

sum = 0;
for i=1:size(N,1)
    xPh = N(i,:)*wcoord;
    fPh(i,1) = subst(f, x, xPh);
    sum = sum + w*N(i,:).'*fPh(i);
end

Fe = he/2*sum;

% xPh1 = N(1,:)*wcoord;
% f1 = subst(f, x, xPh1);
% xPh2 = N(1,:)*wcoord;
% f2 = subst(f, x, xPh2);

% Fe = he/2*(w*N(1,:).'*fPh(1)+w*N(2,:).'*fPh(2));
end

