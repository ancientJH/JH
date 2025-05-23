// Gmsh project created on Fri Sep 27 21:54:18 2024
ell = 0.01;
lc1 = ell*5;
lc2 = ell/2;
dh = 6e-3;
dl = 3e-2;
xlength = 1;
ylength = 1;
Point(1) = {0, 0, 0};
Point(2) = {xlength, 0, 0};
Point(3) = {xlength, ylength, 0};
Point(4) = {0, ylength, 0};
Point(5) = {0, ylength/2-dh/2, 0};
Point(6) = {0, ylength/2+dh/2, 0};
Point(7) = {xlength/2-dl, ylength/2-dh/2, 0};
Point(8) = {xlength/2-dl, ylength/2+dh/2, 0};
Point(9) = {xlength/2, ylength/2, 0};
Point(10) = {xlength, ylength/2, 0};
//+
Line(1) = {1, 2};
Line(2) = {2, 10};
Line(3) = {10, 3};
Line(4) = {3, 4};
Line(5) = {4, 6};
Line(6) = {6, 8};
Line(7) = {8, 9};
Line(8) = {9, 7};
Line(9) = {7, 5};
Line(10) = {5, 1};
// Line(11) = {9, 10};
//+
Circle(12) = {9, 10, 2};
//+
Curve Loop(1) = {10, 1, 2, 3, 4, 5, 6, 7, 8, 9};
Plane Surface(1) = {1};

Line{12} In Surface{1};
Physical Curve("top", 14) = {4};
Physical Curve("bottom", 15) = {1};
Physical Curve("left", 16) = {10, 5};
Physical Curve("right", 17) = {2, 3};
Physical Surface("Surf", 13) = {1};
MeshSize{1:12} = lc1;
MeshSize{2,7,9} = lc2;
Mesh 2;


