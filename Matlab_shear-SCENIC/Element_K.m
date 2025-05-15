function k_e = Element_K(Coord,localCoord, localSol, ElementType, ...
    constitutive, body_force, traction, current_load_step, nQuad, ielem, Psi_plus_rec, Psi_plus_old)%把一个单元的参数传进来
if (nargin < 15)
    switch ElementType
        case 'P12D'
            nQuad = 3; nQuadBdry = 2;%有三个边，每个边有两个点
        case 'Q12D'
            nQuad = 4; nQuadBdry = 2;
    end
end
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

k_e_uu = zeros(length(nDim*nNodesElement));
k_e_dd = zeros(length(nNodesElement));
k_e_ud = zeros(length(nDim*nNodesElement), length(nNodesElement));%还有ud?
k_e_du = zeros(length(nNodesElement), length(nDim*nNodesElement));%还有ud?

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
    StrainVector = B * localSol(u_DoFs);%应力
    d = localSol(d_DoFs)' * Na(:,iQuad);
    %grad_d = Bd * localSol(d_DoFs);
    [~, D, ~, sigma_plus, gc, ell] = constitutive(StrainVector, d);%出来相场的东西
     %BodyForces = current_load_step * body_force(localCoord*Na(:,iQuad));%应该是要转化为高斯点 修改前 两个方向的bodyforce,施加在了高斯点上
      %BodyForces = current_load_step * body_force(Coord);%修改后
      %Psi_plus_rec(iQuad,ielem) = max(Psi_plus_old(iQuad,ielem),Psi_plus);
      Psi_plus_new =  Psi_plus_rec(iQuad,ielem);
    k_e_uu = k_e_uu + B' * D * B * w(iQuad) * detJ(iQuad);
    k_e_dd = k_e_dd + (Na(:,iQuad) * (2 * Psi_plus_new + gc / ell) * Na(:,iQuad)' ...
        + Bd' * gc * ell * Bd) ...
        * w(iQuad) * detJ(iQuad);%修改前
    k_e_ud = k_e_ud - 2 * B' * (1-d) * sigma_plus * Na(:,iQuad)' ...
        * w(iQuad) * detJ(iQuad);%论文中的Cap  修改前
    if Psi_plus_rec(iQuad,ielem)>Psi_plus_old(iQuad,ielem)
        k_e_du = k_e_du - (2 * B' * (1-d) *  sigma_plus * Na(:,iQuad)')' ...
            * w(iQuad) * detJ(iQuad);
    end

end

k_e = zeros(length(localSol), length(localSol));
k_e(u_DoFs, u_DoFs) = k_e_uu;
k_e(d_DoFs, d_DoFs) = k_e_dd;
k_e(u_DoFs, d_DoFs) = k_e_ud;
k_e(d_DoFs, u_DoFs) = k_e_du;