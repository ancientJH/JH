
ell = 0.1e-3;
lc1 = ell/3;
lc2 = ell;
lc3 = ell*5;
L = 25e-3;
H = 5e-3;
L1 = 0.6*H;
H1 = 0.6*H;
//+
Point(1) = {0, 0, 0};
Point(2) = {L, 0, 0};
Point(3) = {L, H, 0};
Point(4) = {0, H, 0};
Point(5) = {L1, H1, 0};
Point(6) = {L, H1, 0};
Point(7) = {L1, H, 0};
//+
Line(1) = {1, 2};
Line(2) = {2, 6};
Line(3) = {6, 3};
Line(4) = {3, 7};
Line(5) = {7, 4};
Line(6) = {4, 1};
Line(7) = {5, 6};
Line(8) = {5, 7};
//+
Curve Loop(1) = {6, 1, 2, 3, 4, 5};
Plane Surface(1) = {1};
Line{7,8} In Surface{1};
Physical Curve("top", 9) = {4, 5};
Physical Curve("bottom", 10) = {1};
Physical Curve("left", 11) = {6};
Physical Curve("right", 12) = {3, 2};

MeshSize{1:4} = lc1;
MeshSize{5:7} = lc2;
MeshSize{3} = lc3;
Mesh 2;//+

//+

