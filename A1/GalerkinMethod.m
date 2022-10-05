function [u] = GalerkinMethod(N,f,b,g,L,rho)
syms x
N_L = subs(N,L); 
B = diff(N,x);
BtB = B.'*B; 
NtN = N.'*N;
K = int(BtB,0,L) - int(NtN,0,L)*rho;
F = int(N.'*f,0,1) + N_L.'*b;
r = 1; l = 2:length(N);

% Solution
dl = K(l,l)\(F(l)+K(l,r)*g); % Originally: - K(l,r)*-g --> dr = -g (boundary)
u = -g+N(l)*dl; % Same as creating one vector d and doing N*d without the -g term.
end