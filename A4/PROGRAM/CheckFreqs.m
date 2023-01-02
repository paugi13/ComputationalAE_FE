% Code to check natural frequencies
E = 70e9;
L=2;
hy=0.25;
e=0.05;

%Number of data points 
I = 1/12*(hy^4-(hy-2*e)^4); % Beam's Section Inertia
V = L*(hy^2-(hy-2*e)^2);
rho = 2700;
C1 = 3.52;
C2 = 4.694^2;
C3 = 7.885^2;
m = V*rho;
m0 = m/L;

f1 = C1/(2*pi)*sqrt(E*I/(m0*L^4));
f2 = C2/(2*pi)*sqrt(E*I/(m0*L^4));
f3 = C3/(2*pi)*sqrt(E*I/(m0*L^4));


