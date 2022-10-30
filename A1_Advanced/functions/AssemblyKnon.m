function [K] = AssemblyKnon(COOR,CN,d_k,AreaFUN,DerStressFUN)
% COOR: Coordinate matrix
% CN % Element connectivy matrix  ;
nelem = size(CN,1) ; % Number of rows and columns SIZE
nnode = nelem+1 ;
nnodeE = size(CN,2) ;
K =zeros(nnode, nnode) ;
for e=1:nelem  % Loop over number of elements 
    % Element matrix
    NODOSe = CN(e,:);    % Global numbering of nodes of element "e"
    COOR_e = COOR(NODOSe) ;
    he = COOR_e(2)-COOR_e(1) ; % Size finite element
    Be = 1/he*[-1 1];
    % Young Modulus & Area Calculation
    epsilon_e = Be*d_k(NODOSe);
    Ee = DerStressFUN(epsilon_e);
    A1 = AreaFUN(COOR_e(1));
    A2 = AreaFUN(COOR_e(2));
    Ae = (A1+A2)/2;
    % Elemental matrix
    Ke = Ee*Ae/he*[1 -1; -1 1]; 
    % Assembly
    for a = 1:nnodeE
        for b = 1:nnodeE
        % Same working principle as the aerospace structures' one.
            A = CN(e,a);
            B = CN(e,b);
            K(A,B) = K(A,B) + Ke(a,b);
        end
    end
end