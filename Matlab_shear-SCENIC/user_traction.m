function Tbar=user_traction(Coord,localCoord)
p=1;
tol=1e-8;
Tbar = zeros(size(Coord));
if  Coord(1,:) - 0<= tol&&abs(Coord(2,:) - 0) <= tol
    if localCoord(2,1)-0>-tol && localCoord(2,2)-0>-tol && localCoord(2,3)-0>-tol
        Tbar(2,1) = -1 * p;
    else
        Tbar(2,1) = 1 * p;
    end
end
% Tbar(2, abs(Coord(1,:) - 10.0) <= 1e-6) = -1.;
