function ValidationTest(F, React, COOR)
% Function to validate obtained results by checking reations values. Both
% Forces and moments are checked.

% COOR = [x1 y1 z1
%       x2 y2 z2 ...]
% F = [F11 F12 F13 F21 F22 F23 ...]
ToL = 10^-4;
% Total reactions values. Initialization.
nnode = size(COOR, 1);
ndim = size(COOR, 2);
xC = 1:3:nnode*ndim;
yC = 2:3:nnode*ndim;
zC = 3:3:nnode*ndim;
Rx = 0; Fx = 0;
Ry = 0; Fy = 0;
Rz = 0; Fz = 0;
Mx = [0 0 0]; MAx = [0 0 0];
My = [0 0 0]; MAy = [0 0 0];
Mz = [0 0 0]; MAz = [0 0 0];
for i = 1:nnode
    xCoor = [COOR(i,1) 0 0];
    yCoor = [0 COOR(i,2) 0];
    zCoor = [0 0 COOR(i,3)];
    % Total acting forces
    Fx = Fx + F(xC(i));
    Fy = Fy + F(yC(i));
    Fz = Fz + F(zC(i));
    %Total acting moments
    FVectorx = [F(xC(i)) 0 0];
    FVectory = [0 F(yC(i)) 0];
    FVectorz = [0 0 F(zC(i))];
    MAx = MAx + cross(yCoor, FVectorz) + cross(zCoor, FVectory);
    MAy = MAy + cross(xCoor, FVectorz) + cross(zCoor, FVectorx);
    MAz = MAz + cross(xCoor, FVectory) + cross(zCoor, FVectorx);
    % Total reaction forces
    Rx = Rx + React(xC(i));
    Ry = Ry + React(yC(i));
    Rz = Rz + React(zC(i));
    %Total reaction moments
    RVectorx = [React(xC(i)) 0 0];
    RVectory = [0 React(yC(i)) 0];
    RVectorz = [0 0 React(zC(i))];
    Mx = Mx + cross(yCoor, RVectorz) + cross(zCoor, RVectory);
    My = My + cross(xCoor, RVectorz) + cross(zCoor, RVectorx);
    Mz = Mz + cross(xCoor, RVectory) + cross(yCoor, RVectorx);
end

% Check results
if Fx + Rx > ToL || Fy + Ry > ToL || Fz + Rz > ToL
    error('Forces are not compensated');
end

if abs(MAx(1) + Mx(1)) > ToL || abs(MAy(2) + My(2)) > ToL || abs(MAz(3) + ...
        Mz(3)) > ToL
    error('Moments are not compensated');
end
   