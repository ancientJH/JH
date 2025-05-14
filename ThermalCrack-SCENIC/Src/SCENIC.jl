function Masterstep_u!(u, ϕ, θ, K, R_u, R_ϕ, R_θ, cellvalues_u, cellvalues_ϕ, cellvalues_θ, dh, dh_u, dh_ϕ, dh_θ, ch, ch_u, ch_ϕ, ch_θ, material, states, states_old, flag_eta)
    @inbounds begin
        u_free = free_dofs(ch_u)
        ϕ_free = free_dofs(ch_ϕ)
        n_u = ndofs(dh_u)
        n_ϕ = ndofs(dh_ϕ)


        #= q = vcat(u, ϕ, θ)  # 合并位移和相场的自由度
        K = assemble_global_K(q, K, cellvalues_u, cellvalues_ϕ, cellvalues_θ, dh, material, states, states_old)
        apply_zero!(K, r, ch);
        Kuu = @view K[1:n_u, 1:n_u]
        Kuϕ = @view K[1:n_u, n_u+1:n_u+n_ϕ]
        Kϕu = @view K[n_u+1:n_u+n_ϕ, 1:n_u]
        Kϕϕ = @view K[n_u+1:n_u+n_ϕ, n_u+n_u+n_ϕ] =#

        q = vcat(u, ϕ)  # 合并位移和相场的自由度
        K = create_sparsity_pattern(dh)
        K = assemble_global_K_ud(q, K, cellvalues_u, cellvalues_ϕ, cellvalues_θ, dh, material, states, states_old)
        apply_zero!(K, vcat(R_u, R_ϕ), ch);

        Kuu = @view K[1:n_u, 1:n_u]
        Kuϕ = @view K[1:n_u, n_u+1:end]
        Kϕu = @view K[n_u+1:end, 1:n_u]
        Kϕϕ = @view K[n_u+1:end, n_u+1:end]
       #=  ################################################################# Kud
        aa = 0.001
        del_ϕ = aa*ones(n_ϕ);
        states_temp = deepcopy(states); states_old_temp = deepcopy(states_old);
        for j in 1:size(states, 2)
            for i in 1:size(states, 1)
                states[i, j] = HistoryVariable(states[i, j].H, states[i, j].ϕ + aa, states[i, j].θ)
            end
        end
        Ru_esti = R_u + Kuϕ * del_ϕ;
        Ru_new = assemble_global_r(u, cellvalues_u, dh_u, material, states_temp, states_old_temp, "u", false);
        println("error = ", norm(Ru_esti - Ru_new)/norm(Ru_esti), "\n")
        ####################################################################    =#
        K_uu = @view (Kuu[u_free, u_free])
        K_uϕ = @view (Kuϕ[u_free, ϕ_free])
        K_ϕu = @view (Kϕu[ϕ_free, u_free])
        K_ϕϕ = @view (Kϕϕ[ϕ_free, ϕ_free])
        

        # Ru, states = assemble_global_r!(u, ϕ, cellvalues_u, cellvalues_ϕ, dh, dh_u, dh_ϕ, material, states, states_old, "u")
        # Rϕ, _ = assemble_global_r!(u, ϕ, cellvalues_u, cellvalues_ϕ, dh, dh_u, dh_ϕ, material, states, states_old, "ϕ")

        K_uu = sparse(Symmetric(K_uu))
        K_ϕϕ = sparse(Symmetric(K_ϕϕ))
        # K_uu_fact = cholesky(K_uu)
        # K_ϕϕ_fact = cholesky(K_ϕϕ)
        # R_u = Ru[u_free]
        # R_ϕ = Rϕ[ϕ_free]
        eta = 1;
        if flag_eta == 1
            rho = cal_rho1(K_uu, K_ϕϕ, K_uϕ, K_ϕu)
            
            if rho < 1
                eta = 1;
                println("eta2 = ", eta)
            elseif rho > 1
                eta = 0.9/sqrt(rho);
                println("eta2 = ", eta)
            end
        end

        y_active = -(K_uu \ (R_u[u_free] - eta * K_uϕ * (K_ϕϕ \ R_ϕ[ϕ_free])))
        #= rho = cal_rho1(K_uu, K_ϕϕ, K_uϕ, K_ϕu)
        println("rho = ", rho) =#

        Δu_active = gmres_user(K_uu, K_ϕϕ, eta * K_uϕ, eta * K_ϕu, y_active; max_iter=50, tol=1e-6)
        Δu = zeros(ndofs(dh_u))
        #println("Δu_active自由度为",length(Δu_active))
        #println("Δu[u_free]自由度为",length(Δu[u_free]))
        Δu[u_free] .= Δu_active
    end
    return Δu
end

function Masterstep_θ!(u, ϕ, θ, K, R_u, R_ϕ, R_θ, cellvalues_u, cellvalues_ϕ, cellvalues_θ, dh, dh_u, dh_ϕ, dh_θ, ch, ch_u, ch_ϕ, ch_θ, material, states, states_old, flag_eta)
    @inbounds begin
        u_free = free_dofs(ch_u)
        ϕ_free = free_dofs(ch_ϕ)
        θ_free = free_dofs(ch_θ)
        n_u = ndofs(dh_u)
        n_ϕ = ndofs(dh_ϕ)
        n_θ = ndofs(dh_θ)
        q = vcat(u, ϕ, θ)  # 合并位移和相场的自由度
        K = assemble_global_K(q, K, cellvalues_u, cellvalues_ϕ, cellvalues_θ, dh, material, states, states_old)

        Kuu = @view K[1:n_u, 1:n_u]
        Kuϕ = @view K[1:n_u, n_u+1:n_u+n_ϕ]
        Kuθ = @view K[1:n_u, n_u+n_ϕ+1:end]
        Kϕu = @view K[n_u+1:n_u+n_ϕ, 1:n_u]
        Kϕϕ = @view K[n_u+1:end, n_u+1:end]
        Kϕθ = @view K[n_u+1:end, n_u+n_ϕ+1:end]
        Kθu = @view K[n_u+n_ϕ+1:end, 1:n_u]
        Kθϕ = @view K[n_u+n_ϕ+1:end, n_u+1:n_u+n_ϕ]
        Kθθ = @view K[n_u+n_ϕ+1:end, n_u+n_ϕ+1:end]
        K_uu = @view (Kuu[u_free, u_free])
        K_uϕ = @view (Kuϕ[u_free, ϕ_free])
        K_uθ = @view (Kuθ[u_free, θ_free])
        K_ϕu = @view (Kϕu[ϕ_free, u_free])
        K_ϕϕ = @view (Kϕϕ[ϕ_free, ϕ_free])
        K_ϕθ = @view (Kϕθ[ϕ_free, θ_free])
        K_θu = @view (Kθu[θ_free, u_free])
        K_θϕ = @view (Kθϕ[θ_free, ϕ_free])
        K_θθ = @view (Kθθ[θ_free, θ_free])

        K_uu = sparse(Symmetric(K_uu))
        K_ϕϕ = sparse(Symmetric(K_ϕϕ))
        K_θθ = sparse(Symmetric(K_θθ))

        eta = 1;
        if flag_eta == 1
            rho = cal_rho1(
                            K_θθ,                                                               # 对应原 K_uu
                            vcat(hcat(K_uu, K_uϕ), hcat(K_ϕu, K_ϕϕ)),                            # 对应原 K_ϕϕ
                            hcat(K_θu, K_θϕ),                                                     # 对应原 K_uϕ
                            vcat(K_uθ, K_ϕθ)                                                      # 对应原 K_ϕu
                        )
            if rho < 1
                eta = 1;
                println("eta1 = ", eta)
            elseif rho > 1
                eta = 0.9/sqrt(rho);
                println("eta1 = ", eta)
            end
        end

        y_active = -(K_θθ \ (R_θ[θ_free] - eta * hcat(K_θu, K_θϕ) * (vcat(hcat(K_uu, K_uϕ), hcat(K_ϕu, K_ϕϕ)) \ vcat(R_u[u_free], R_ϕ[ϕ_free]))))

        Δθ_active = gmres_user(K_θθ, vcat(hcat(K_uu, K_uϕ), hcat(K_ϕu, K_ϕϕ)), eta * hcat(K_θu, K_θϕ), eta * vcat(K_uθ, K_ϕθ), y_active; max_iter=50, tol=1e-6)
        Δθ = zeros(ndofs(dh_θ))
        Δθ[θ_free] .= Δθ_active
    end
    return Δθ
end

function gmres_user(K11, K22, K12, K21, b; max_iter=50, tol=1e-6)
    n = length(b)
    x0 = zeros(n)
    r0 = b - closure_M(K11, K22, K12, K21, x0)
    β = norm(r0)
    

    V = zeros(n, max_iter + 1)
    V[:, 1] = r0 / β
    ξ = zeros(max_iter + 1)
    ξ[1] = β
    H = zeros(max_iter + 1, max_iter)
    c = zeros(max_iter)
    s = zeros(max_iter)

    j = 0
    for j in 1:max_iter
        # Arnoldi 过程
        w = closure_M(K11, K22, K12, K21, V[:, j])
        for i in 1:j
            H[i, j] = dot(V[:, i], w)
            w -= H[i, j] * V[:, i]
        end
        for i in 1:j
            h_ij = dot(V[:, i], w)
            H[i, j] += h_ij
            w -= h_ij * V[:, i]
        end
        H[j+1, j] = norm(w)
        if H[j+1, j] < eps() * β
            break
        end
        V[:, j+1] = w / H[j+1, j]

        # Givens 旋转
        for i in 1:j-1
            temp = c[i] * H[i, j] + s[i] * H[i+1, j]
            H[i+1, j] = -s[i] * H[i, j] + c[i] * H[i+1, j]
            H[i, j] = temp
        end
        rho = sqrt(H[j, j]^2 + H[j+1, j]^2)
        if rho != 0
            c[j] = H[j, j] / rho
            s[j] = H[j+1, j] / rho
            H[j, j] = rho
            H[j+1, j] = 0.0
        else
            c[j] = 1.0
            s[j] = 0.0
        end

        # 更新 ξ
        ξ[j+1] = -s[j] * ξ[j]
        ξ[j] = c[j] * ξ[j]

        relres = abs(ξ[j+1]) / β
        if relres < tol 
            y = UpperTriangular(H[1:j, 1:j]) \ ξ[1:j]
            x = x0 + V[:, 1:j] * y
            # println("第 $j 次迭代收敛，相对残差为 $relres")
            return x
        end
    end

    # 如果达到最大迭代次数或需要最终解
    y = UpperTriangular(H[1:j, 1:j]) \ ξ[1:j]
    x = x0 + V[:, 1:j] * y
    relres = abs(ξ[j+1]) / β
    error("GMRES 未收敛，迭代 $j 次，相对残差为 $relres")
    return x
end

function cal_rho1(K_uu, K_ϕϕ, K_uϕ, K_ϕu)
    n = size(K_uu, 1)  # K_uu 的行数
    m = size(K_ϕϕ, 1)  # K_ϕϕ 的行数
    A_top = hcat(zeros(n, n), -K_uϕ)        # 第一行块: [0 | -K_uϕ]
    A_bottom = hcat(zeros(m, n), zeros(m, m)) # 第二行块: [0 | 0]
    A = vcat(A_top, A_bottom)

    # 构造完整矩阵 B（右侧矩阵）
    B_top = hcat(K_uu, zeros(n, m))          # 第一行块: [K_uu | 0]
    B_bottom = hcat(K_ϕu, K_ϕϕ)              # 第二行块: [K_vu | K_vv]
    B = vcat(B_top, B_bottom)
    rho = eigs(sparse(A), sparse(B); nev=1, which=:LM, tol=1e-2)[1]  # 只提取模最大的特征值
    rho = abs.(rho)
    rho = Float64(rho[])
    return rho
end

function cal_rho(K_uu, K_ϕϕ, K_uϕ, K_ϕu, R_u) # 老办法计算rho
    max_iter = 20  # 最大迭代次数
    tol = 1e-2    # 收敛容差
    R_u_old = R_u
    rho = 0.0     # 初始化 rho
    prev_rho = 0.0  # 上一次的 rho 值
    for i ∈ 1:max_iter
        R_u_new = K_uu \ (K_uϕ * (K_ϕϕ \ (K_ϕu * R_u_old)))
        rho = norm(R_u_new) / norm(R_u_old)
        if abs(rho - prev_rho) < tol
            println("cal_rho: Converged at iteration $i with rho_diff = $(abs(rho - prev_rho))")
            return rho
        end
        prev_rho = rho
        R_u_old = R_u_new
    end
    println("cal_rho: Reached max_iter without convergence")
    return rho
end

function closure_M(K_uu, K_ϕϕ, K_uϕ, K_ϕu, x)
    result = x - K_uu\(K_uϕ*(K_ϕϕ\(K_ϕu*x)));
    return result
end

function updatePsi(cellvalues_u, u, dh_u, material, states, states_old)
    for (i, cell) in enumerate(CellIterator(dh_u))
        reinit!(cellvalues_u, cell)
        eldofs_u = celldofs(cell)         # 针对dh_u，编号就是u场的本地编号
        ue = u[eldofs_u]
        state = @view states[:, i]
        state_old = @view states_old[:, i]
        states[:, i] = assemble_element_psi!(cellvalues_u, ue, material, state, state_old)
    end
    return states
end

function assemble_element_psi!(cellvalues_u,ue::Vector, material::Material, state, state_old)
    for q_point in 1:getnquadpoints(cellvalues_u)
        ϕ = state[q_point].ϕ;
        θ = state[q_point].θ;
        εθ = material.α*(θ - material.θ₀);
        if material.dim == 2
            ε_PlaneStrain = function_symmetric_gradient(cellvalues_u, q_point, ue) - εθ*one(SymmetricTensor{2, 2});
            ∇ˢu = SymmetricTensor{2,3,Float64}((ε_PlaneStrain[1, 1], ε_PlaneStrain[1, 2], 0.0, ε_PlaneStrain[2, 2], 0.0, 0.0));
        elseif material.dim == 3
            ∇ˢu = function_symmetric_gradient(cellvalues_u, q_point, ue) - εθ*one(SymmetricTensor{2, 3});
        else
            error("Invalid element dimension");
        end
        
        gd, _, _ = Degradation(ϕ, material);
        Ψ⁺, _, σ, _ = Constitutive(∇ˢu, material, gd);
        H = max(Ψ⁺, state_old[q_point].H);
        state[q_point] = HistoryVariable(H, ϕ, θ);
    end
    return state
end

function updateT(cellvalues_θ, θ, dh_θ, material, states, states_old)
    for (i, cell) in enumerate(CellIterator(dh_θ))
        reinit!(cellvalues_θ, cell)
        eldofs_θ = celldofs(cell)         
        θe = θ[eldofs_θ]
        state = @view states[:, i]
        state_old = @view states_old[:, i]
        states[:, i] = assemble_element_T!(cellvalues_θ, θe, material, state, state_old)
    end
    return states
end

function assemble_element_T!(cellvalues_θ,θe::Vector, material::Material, state, state_old)
    for q_point in 1:getnquadpoints(cellvalues_θ)
        H = state[q_point].H;
        ϕ = state[q_point].ϕ;
        θ = function_value(cellvalues_θ, q_point, θe);
        state[q_point] = HistoryVariable(H, ϕ, θ);
    end
    return state
end
