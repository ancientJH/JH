function R = Assemble_Rd(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
    constitutive, body_force, traction, current_load_step, Psi_plus_rec, Psi_plus_old)
nNodesElement = size(IEN, 1);%每个单元有几个节点
nElements = size(IEN, 2);
nDim = size(Coord,1);
nNodes = size(Coord, 2);
nEquations_u = length(Sol_u);
nEquations_d = length(Sol_d);
nDoF = nDim+1;
R = zeros(size(Sol_d));
switch elementType
    case 'P12D'
        nQuad = 3; nQuadBdry = 2;   %Nquad is the nodes'number of one element,nQuadRdry is the nodes'number of one side
    case 'Q12D'
        nQuad = 4; nQuadBdry = 2;
end
local_d_indices = nDoF : nDoF : (nDoF * nNodesElement);
local_u_indices = setdiff(1:(nDoF * nNodesElement), local_d_indices);%C = setdiff(A,B) 返回 A 中存在但 B 中不存在的数据，...
...不包含重复项。C 是有序的。 对于一个三角形单元，每个节点有两个位移量，一个相场量，位移量占前两个，相场量占最后一个
for ielem = 1:nElements
    localCoord = Coord(:, IEN(:, ielem));%第i个单元的三个节点的坐标
    localSol_u = Sol_u(LM_u(:, ielem));
    localSol_d = Sol_d(LM_d(:, ielem));
    localSol = zeros(length(localSol_u) + length(localSol_d));
    localSol(local_u_indices) = localSol_u;
    localSol(local_d_indices) = localSol_d;
    r_e = Element_Rd(localCoord, localSol, elementType, ...
    constitutive, body_force, current_load_step, nQuad, ielem, Psi_plus_rec, Psi_plus_old);
    R(LM_d(:, ielem)) = R(LM_d(:, ielem)) + r_e;
end


