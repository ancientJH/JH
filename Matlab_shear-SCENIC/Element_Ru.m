function [r_e, Psi_plus_rec] = Element_Ru(localCoord, localSol, ElementType, ...
    constitutive, body_force, current_load_step, nQuad, ielem, Psi_plus_rec, Psi_plus_old)%把一个单元的参数传进来
if (nargin < 15)
    switch ElementType
        case 'P12D'
            nQuad = 3; nQuadBdry = 2;%有三个边，每个边有两个点
        case 'Q12D'
            nQuad = 4; nQuadBdry = 2;
    end
end
% global Psi_plus_rec Psi_plus_old
[xi, w] = GetQuadratureRule(ElementType, nQuad);
nDim = size(localCoord, 1);%二维问题
nDoF = nDim + 1;%每个节点有三个量：两个位移+1个d
nNodesElement = size(localCoord, 2);

d_DoFs = nDoF * (1:nNodesElement)';
u_DoFs = true(length(localSol), 1);
u_DoFs(d_DoFs) = false;

switch nDim
    case 1
        [detJ, Na, dNa_dx] = QuadShape(ElementType, localCoord, xi);
    case 2
        [detJ, Na, dNa_dx, dNa_dy] = QuadShape(ElementType, localCoord, xi);
    case 3
        [detJ, Na, dNa_dx, dNa_dy, dNa_dz] = QuadShape(ElementType, localCoord, xi);
end
r_e_u = zeros(nDim*nNodesElement, 1);

zero_filler = zeros(1, length(dNa_dx));
for iQuad = 1:nQuad
    switch nDim
        case 1
            B = dNa_dx(:,iQuad);
            Bd = B;
        case 2
            B = [reshape([dNa_dx(:,iQuad)'; zero_filler], 1, [])
                reshape([zero_filler; dNa_dy(:,iQuad)'], 1, [])
                reshape([dNa_dy(:,iQuad)'; dNa_dx(:,iQuad)'], 1, [])];%为什么要这样组合？
            Bd = [dNa_dx(:,iQuad)'; dNa_dy(:,iQuad)'];%为什么要这样组合？
        case 3
            B = [reshape([dNa_dx(:,iQuad)'; zero_filler; zero_filler], 1, [])
                reshape([zero_filler; dNa_dy(:,iQuad)'; zero_filler], 1, [])
                reshape([zero_filler; zero_filler; dNa_dz(:,iQuad)'], 1, [])
                reshape([zero_filler; dNa_dz(:,iQuad)'; dNa_dy(:,iQuad)'], 1, [])
                reshape([dNa_dz(:,iQuad)'; zero_filler; dNa_dx(:,iQuad)'], 1, [])
                reshape([dNa_dy(:,iQuad)'; dNa_dx(:,iQuad)'; zero_filler], 1, [])];
            Bd = [dNa_dx(:,iQuad)'; dNa_dy(:,iQuad)'; dNa_dz(:,iQuad)'];
    end
    d = localSol(d_DoFs)' * Na(:,iQuad);
    StrainVector = B * localSol(u_DoFs);%应力
    [StressVector, ~, Psi_plus, ~, ~, ~] = constitutive(StrainVector, d);%出来相场的东西
    Psi_plus_rec(iQuad,ielem) = max(Psi_plus_old(iQuad,ielem),Psi_plus);
    BodyForces = current_load_step * body_force(localCoord*Na(:,iQuad));
    r_e_u = r_e_u + (B' * StressVector - reshape(BodyForces * Na(:,iQuad)', [], 1)) ...
        * w(iQuad) * detJ(iQuad);
end

switch ElementType
    case 'P12D'
        nFaces = 3;
        FaceNodes = [1, 2, 3
            2, 3, 1];
        FaceElementType = 'P11D';
    case 'Q12D'
        nFaces = 4;
        FaceNodes = [1, 2, 3, 4
            2, 3, 4, 1];
        FaceElementType = 'P11D';
    case 'P13D'
        nFaces = 4;
        FaceNodes = [1, 1, 1, 2
            3, 2, 4, 3
            2, 4, 3, 4];
end

% Compute the contribution of the traction to the residual.
[xi, w] = GetQuadratureRule(FaceElementType, nQuadBdry);
for iFace = 1:nFaces
    ThisFaceNodes = FaceNodes(:, iFace);
    FaceCoord = localCoord(:, ThisFaceNodes);
    element_uDoFs = reshape(repmat(ThisFaceNodes'-1, nDim, 1) * nDim + repmat((1:nDim)', 1, length(ThisFaceNodes)), [], 1);
    [detJ, Na] = FaceQuadShape(FaceElementType, FaceCoord, xi);
    for iQuadBdry = 1:nQuadBdry
%         TractionThisPoint = current_load_step * traction(FaceCoord *
%         Na(:, iQuadBdry));修改前
%TractionThisPoint = current_load_step * traction(FaceCoord * Na(:, iQuadBdry),localCoord);
 TractionThisPoint = [0;0];
        r_e_u(element_uDoFs) = r_e_u(element_uDoFs) - reshape(TractionThisPoint * Na(:, iQuadBdry)', [], 1) ...
            * w(iQuadBdry) * detJ(iQuadBdry);
    end
end

r_e = r_e_u;

