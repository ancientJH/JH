function Kuu = Assemble_Ku(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
    constitutive, body_force, traction, current_load_step, Psi_plus_rec, Psi_plus_old)
    
    % 基本参数初始化
    nNodesElement = size(IEN, 1);
    nElements = size(IEN, 2);
    nDim = size(Coord, 1);
    nDoF = nDim + 1;
    
    % 自由度索引定义
    local_d_indices = nDoF : nDoF : (nDoF * nNodesElement);
    local_u_indices = setdiff(1:(nDoF * nNodesElement), local_d_indices);
    %n_u = length(local_u_indices);
    %n_d = length(local_d_indices);
    
    % 预计算全局矩阵非零元素估计
    avg_nnz_per_element = nNodesElement*nDoF;
    total_nnz = ceil(1.5 * avg_nnz_per_element^2 * nElements);
    
    % 初始化三元组存储数组
    [i_uu, j_uu, v_uu] = deal(zeros(total_nnz, 1));
    
    ptr = ones(1, 1);
    
    % 积分规则选择
    switch elementType
        case 'P12D'
            nQuad = 3; nQuadBdry = 2;
        case 'Q12D'
            nQuad = 4; nQuadBdry = 2;
        otherwise
            error('Unsupported element type');
    end
    
    % 主组装循环
    for ielem = 1:nElements
        localCoord = Coord(:, IEN(:, ielem));
        u_nodes = LM_u(:, ielem);
        d_nodes = LM_d(:, ielem);
        
        localSol = zeros(nDoF * nNodesElement, 1);
        localSol(local_u_indices) = Sol_u(u_nodes);
        localSol(local_d_indices) = Sol_d(d_nodes);
        
        k_uu = Element_Ku(Coord, localCoord, localSol, elementType, constitutive, ...
            body_force, traction, current_load_step, nQuad, ielem, Psi_plus_rec, Psi_plus_old);
        
        
        
        
        % 修正索引生成部分
        % Kuu
        [I, J] = ndgrid(u_nodes, u_nodes);
        idx = ptr(1):(ptr(1)+numel(k_uu)-1);
        i_uu(idx) = I(:);
        j_uu(idx) = J(:);
        v_uu(idx) = k_uu(:);
        ptr(1) = ptr(1) + numel(k_uu);

    end
    Kuu = sparse(i_uu(1:ptr(1)-1), j_uu(1:ptr(1)-1), v_uu(1:ptr(1)-1), length(Sol_u), length(Sol_u));
end