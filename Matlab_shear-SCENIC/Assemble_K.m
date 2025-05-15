function [Kuu, Kdd, Kud, Kdu] = Assemble_K(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
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
    [i_dd, j_dd, v_dd] = deal(zeros(total_nnz, 1));
    [i_ud, j_ud, v_ud] = deal(zeros(total_nnz, 1));
    [i_du, j_du, v_du] = deal(zeros(total_nnz, 1));
    ptr = ones(4, 1);
    
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
        
        k_e = Element_K(Coord, localCoord, localSol, elementType, constitutive, ...
            body_force, traction, current_load_step, nQuad, ielem, Psi_plus_rec, Psi_plus_old);
        
        k_uu = k_e(local_u_indices, local_u_indices);
        k_dd = k_e(local_d_indices, local_d_indices);
        k_ud = k_e(local_u_indices, local_d_indices);
        k_du = k_e(local_d_indices, local_u_indices);
        
        % 修正索引生成部分
        % Kuu
        [I, J] = ndgrid(u_nodes, u_nodes);
        idx = ptr(1):(ptr(1)+numel(k_uu)-1);
        i_uu(idx) = I(:);
        j_uu(idx) = J(:);
        v_uu(idx) = k_uu(:);
        ptr(1) = ptr(1) + numel(k_uu);
        
        % Kdd
        [I, J] = ndgrid(d_nodes, d_nodes);
        idx = ptr(2):(ptr(2)+numel(k_dd)-1);
        i_dd(idx) = I(:);
        j_dd(idx) = J(:);
        v_dd(idx) = k_dd(:);
        ptr(2) = ptr(2) + numel(k_dd);
        
        % Kud (关键修正)
        [I, J] = ndgrid(u_nodes, d_nodes); % 行是u_nodes，列是d_nodes
        idx = ptr(3):(ptr(3)+numel(k_ud)-1);
        i_ud(idx) = I(:);
        j_ud(idx) = J(:);
        v_ud(idx) = k_ud(:);
        ptr(3) = ptr(3) + numel(k_ud);
        
        % Kdu
        [I, J] = ndgrid(d_nodes, u_nodes); % 行是d_nodes，列是u_nodes
        idx = ptr(4):(ptr(4)+numel(k_du)-1);
        i_du(idx) = I(:);
        j_du(idx) = J(:);
        v_du(idx) = k_du(:);
        ptr(4) = ptr(4) + numel(k_du);
    end
    
    % 截断多余空间并构建稀疏矩阵
    Kuu = sparse(i_uu(1:ptr(1)-1), j_uu(1:ptr(1)-1), v_uu(1:ptr(1)-1), length(Sol_u), length(Sol_u));
    Kdd = sparse(i_dd(1:ptr(2)-1), j_dd(1:ptr(2)-1), v_dd(1:ptr(2)-1), length(Sol_d), length(Sol_d));
    Kud = sparse(i_ud(1:ptr(3)-1), j_ud(1:ptr(3)-1), v_ud(1:ptr(3)-1), length(Sol_u), length(Sol_d));
    Kdu = sparse(i_du(1:ptr(4)-1), j_du(1:ptr(4)-1), v_du(1:ptr(4)-1), length(Sol_d), length(Sol_u));
end