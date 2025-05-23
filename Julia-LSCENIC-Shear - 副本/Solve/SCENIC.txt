function Masterstep_u!(u, ϕ, K_uu, K_ϕϕ, R_u, R_ϕ, cellvalues_u, cellvalues_ϕ, dh, dh_u, dh_ϕ, ch, ch_u, ch_ϕ, material, states, states_old, flag_eta)
    @inbounds begin
        u_free = free_dofs(ch_u)
        ϕ_free = free_dofs(ch_ϕ)
        K_uu = assemble_global_Ku(u, K_uu, cellvalues_u, dh_u, material, states, states_old)
        K_ϕϕ = assemble_global_Kϕ(ϕ, K_ϕϕ, cellvalues_ϕ, dh_ϕ, material, states, states_old)
        K_uu_active = @view K_uu[u_free, u_free];
        K_ϕϕ_active = @view K_ϕϕ[ϕ_free, ϕ_free];
        K_uu_active = sparse(Symmetric(K_uu_active))
        K_ϕϕ_active = sparse(Symmetric(K_ϕϕ_active))
        
        #= ϕ_temp = copy(ϕ); ϕ_temp[ϕ_free] -= K_ϕϕ_active \ R_ϕ;
        r_u_new = assemble_global_ru_withϕ(u, ϕ_temp, cellvalues_u, cellvalues_ϕ, dh, dh_u, dh_ϕ, material, states, states_old, "u", false);
        y_active = -K_uu_active \ (2 * R_u - r_u_new[u_free]);  =#
        y_active = K_uu_active \ (-R_u); 
        Δu_active = gmres_user(K_uu_active, K_ϕϕ_active, R_u, R_ϕ, y_active, u, ϕ, cellvalues_u, cellvalues_ϕ, dh, dh_u, dh_ϕ, material, states, states_old, u_free, ϕ_free, flag_eta; max_iter=50, tol=1e-6)
        Δu = zeros(ndofs(dh_u))
        Δu[u_free] .= Δu_active   
    end
    return Δu
end


function gmres_user(K11, K22, r_u_old, r_ϕ_old, b::Vector{Float64}, u, ϕ, cellvalues_u, cellvalues_ϕ, dh, dh_u, dh_ϕ, Mat, states, states_old, u_free, ϕ_free, flag_eta; max_iter::Int=50, tol::Float64=1e-6)::Vector{Float64}
    @inbounds begin    
        n = length(b)
        x0 = zeros(Float64, n) # GMRES 内部的初始迭代猜测（通常为0）

        r0_vec = b # 初始残差 r_0 = b (假设 A*x0 = 0, 因 x0=0)
        
        norm_b_rhs = norm(b)

        β = norm(r0_vec)

        V_basis = zeros(Float64, n, max_iter + 1)
        H_hessenberg = zeros(Float64, max_iter + 1, max_iter)
        
        cs_givens = zeros(Float64, max_iter)
        sn_givens = zeros(Float64, max_iter)
        
        g_rhs_rotated = zeros(Float64, max_iter + 1)
        g_rhs_rotated[1] = β

        V_basis[:, 1] .= r0_vec ./ β
        u_temp = zeros(Float64, ndofs(dh_u))
        ϕ_temp = zeros(Float64, ndofs(dh_ϕ))
        k_iters = 0
        @inbounds for k in 1:max_iter
            k_iters = k

            current_V_k_view = @view V_basis[:, k]
            V_kplus1_view = @view V_basis[:, k+1]

            # Arnoldi 迭代：计算 A * v_k
            # 使用您指定的 leftterm 参数顺序
            copyto!(u_temp, u); copyto!(ϕ_temp, ϕ);
            # u_temp .= u; ϕ_temp .= ϕ;
            V_kplus1_view .= leftterm(
                u, ϕ, cellvalues_u, cellvalues_ϕ, dh, dh_u, dh_ϕ, Mat, states, states_old,
                current_V_k_view, # 这是 input_vector (即 x0 在您签名中的位置)
                r_u_old, r_ϕ_old, K11, K22, u_temp, ϕ_temp,
                u_free, ϕ_free, flag_eta
            )

            # Gram-Schmidt 正交化
            @inbounds for j_orth in 1:k
                h_coeff = dot((@view V_basis[:, j_orth]), V_kplus1_view)
                H_hessenberg[j_orth, k] = h_coeff
                axpy!(-h_coeff, (@view V_basis[:, j_orth]), V_kplus1_view)
            end

            h_kplus1_k = norm(V_kplus1_view)
            H_hessenberg[k+1, k] = h_kplus1_k

            if h_kplus1_k < eps(Float64) 
                break
            end
            V_kplus1_view ./= h_kplus1_k

            # Givens 旋转
            @inbounds for j_rot in 1:(k-1)
                h_upper = H_hessenberg[j_rot, k]; h_lower = H_hessenberg[j_rot+1, k]
                cs_val = cs_givens[j_rot];       sn_val = sn_givens[j_rot]
                H_hessenberg[j_rot, k]   = cs_val * h_upper + sn_val * h_lower
                H_hessenberg[j_rot+1, k] = -sn_val * h_upper + cs_val * h_lower
            end

            rho_givens = hypot(H_hessenberg[k,k], H_hessenberg[k+1,k])
            local c_k_rot, s_k_rot
            if rho_givens > eps(Float64)
                c_k_rot = H_hessenberg[k,k] / rho_givens
                s_k_rot = H_hessenberg[k+1,k] / rho_givens
                H_hessenberg[k,k] = rho_givens
                H_hessenberg[k+1,k] = 0.0
            else
                c_k_rot = 1.0; s_k_rot = 0.0
            end
            cs_givens[k] = c_k_rot; sn_givens[k] = s_k_rot

            g_rhs_rotated[k+1] = -s_k_rot * g_rhs_rotated[k]
            g_rhs_rotated[k]   =  c_k_rot * g_rhs_rotated[k]
            
            if abs(g_rhs_rotated[k+1]) / norm_b_rhs < tol
                println("GMRES converged in $(k) iterations, the residual is $(abs(g_rhs_rotated[k+1]) / norm_b_rhs).")
                break
            end
        end

        y_solution = UpperTriangular(@view H_hessenberg[1:k_iters, 1:k_iters]) \ (@view g_rhs_rotated[1:k_iters])
        # 最终解 x = x0 + V_m * y_m。因为 x0 是零向量，所以 x = V_m * y_m。
        final_solution = (@view V_basis[:, 1:k_iters]) * y_solution 
    end
    return final_solution
end

function leftterm(u, ϕ, cellvalues_u, cellvalues_ϕ, dh, dh_u, dh_ϕ, Mat, states, states_old, du_active, r_u_old, r_ϕ_old, K_uu, K_ϕϕ, u_temp, ϕ_temp, u_free, ϕ_free, flag_eta) # 传入的应该都是active部分
    e = sqrt(eps(Float64)) * (1.0 + norm(u[u_free])); du_active_scaled = e * du_active; # Renamed to avoid confusion with the original du_active argument
    # u_free = Ferrite.free_dofs(ch_u); # Now passed as argument
    # ϕ_free = Ferrite.free_dofs(ch_ϕ); # Now passed as argument
    u_temp[u_free] .+= du_active_scaled; # Potentially more efficient if du_active_scaled is a vector of correct size
    r_ϕ = assemble_global_rϕ_withu(u_temp, ϕ, cellvalues_u, cellvalues_ϕ, dh, dh_u, dh_ϕ, Mat, states, states_old, "ϕ", true);
    dR_ϕ = r_ϕ[ϕ_free] - r_ϕ_old;
    dϕ_active = K_ϕϕ \ (-dR_ϕ);
    ϕ_temp[ϕ_free] += dϕ_active;

    r_u = assemble_global_ru_withϕ(u, ϕ_temp, cellvalues_u, cellvalues_ϕ, dh, dh_u, dh_ϕ, Mat, states, states_old, "u", true);
    dR_u = r_u[u_free] - r_u_old;
    ddu_active = K_uu \ (-dR_u); 
    eta = 1.0
    if flag_eta == true
        rho = norm(ddu_active)/norm(du_active_scaled) # Use scaled version for norm calculation consistency
        if rho > 1
        eta = 0.9/sqrt(rho)
        end
    end
    return (du_active_scaled - eta^2 * ddu_active)./e # Return based on scaled version
end


function cal_rho1(K_uu, K_ϕϕ, K_uϕ, K_ϕu)
    n = size(K_uu, 1) 
    m = size(K_ϕϕ, 1) 
    A = [spzeros(Float64, n, n)  -K_uϕ;
                spzeros(Float64, m, n)  spzeros(Float64, m, m)]
    B = [K_uu                     spzeros(Float64, n, m);
                K_ϕu                     K_ϕϕ]
    rho = eigs(sparse(A), sparse(B); nev=1, which=:LM, tol=1e-2)[1] 
    rho = abs(rho[1]) 
    return Float64(rho)
end