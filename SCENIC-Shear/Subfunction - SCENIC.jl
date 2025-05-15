#定义材料参数相关的结构体
struct Material{T}
    E::T  # Young's modulus
    G::T  # Shear modulus = μ (Lame constant)
    K::T  # Bulk modulus
    λ::T  # Lame constant
    ν::T  # Poisson's ratio
    Gc::T # Fracture Toughness
    σc::T # Strength
    ℓ::T  # Phase field length scale
    s::T  # parameter for Hughes degradation 
    a₁::T # parameter for Wu degradation 
    flag::StrainDecomp
    flagD::DegradType
    dim::Int64
end
function Material(E, ν, Gc, σc, ℓ, s, flag, flagD, dim)
    G = E / 2(1 + ν);
    K = E / 3(1 - 2ν);
    λ = K - 2G / 3;
    a₁ = 27E*Gc/(128σc^2*ℓ₀);
    return Material(E, G, K, λ, ν, Gc, σc, ℓ, s, a₁, flag, flagD, dim)
end

#定义历史变量结构体
struct HistoryVariable{T}
	H::T # History variable
	ϕ::T #phase field variable from last increment
end
function HistoryVariable()
	return HistoryVariable(0.0, 0.0)
end
#定义求解参数结构体
mutable struct SolverState{T,F}
    loadsteps::Vector{F}
    nitr_inner::T
    nitr_outer::T
    TOL_u::F
    TOL_ϕ::F
end
function create_solver_state(loadsteps, nitr_inner, nitr_outer, TOL_u, TOL_ϕ)
    return SolverState(loadsteps, nitr_inner, nitr_outer, TOL_u, TOL_ϕ)
end
#定义输出参数结构体
mutable struct OutputVariables{T}
    plotframe::T
    totalIterations_outer::T
    totalIterations_ϕ::T
    totalIterations_u::T
    plotFrequency::T
    historyFrequency::T
    a0::Float64
    CrackDir::T
    OutputSet::String
end
#输出参数结构体的初始化
function OutputVariables(field_frequency, history_frequency, a0, CrackDir, outputset)
    return OutputVariables(0, 0, 0, 0, field_frequency, history_frequency, a0, CrackDir, outputset)
end

function assemble_element_u!(Ke, Re, cellvalues_u, ue, material, state, state_old)
    nbase_u = getnbasefunctions(cellvalues_u);
    for q_point in 1:getnquadpoints(cellvalues_u)
        dΩᵤ = getdetJdV(cellvalues_u, q_point)
        if material.dim == 2
            ε_PlaneStrain = function_symmetric_gradient(cellvalues_u, q_point, ue)
            ∇ˢu = SymmetricTensor{2,3,Float64}((ε_PlaneStrain[1, 1], ε_PlaneStrain[1, 2], 0.0,ε_PlaneStrain[2, 2], 0.0, 0.0))
        elseif material.dim == 3
            ∇ˢu = function_symmetric_gradient(cellvalues_u, q_point, ue)
        else
            error("Invalid element dimension")
        end

        ϕ = state[q_point].ϕ
        gd, _, _ = Degradation(ϕ, material);
        Ψ⁺, Dᵉ, σ = Constitutive(∇ˢu, material, gd);
        H = max(Ψ⁺, state_old[q_point].H);
        state[q_point] = HistoryVariable(H, ϕ);

        for i in 1:nbase_u
            if material.dim == 2
                δε_2d = shape_symmetric_gradient(cellvalues_u, q_point, i);
                δε = SymmetricTensor{2,3,Float64}((δε_2d[1, 1], δε_2d[1, 2], 0.0,δε_2d[2, 2], 0.0, 0.0));
            else
                δε = shape_symmetric_gradient(cellvalues_u, q_point, i);
            end
            for j in 1:i
                if material.dim == 2
                    ε_2d = shape_symmetric_gradient(cellvalues_u, q_point, j);
                    ε̄u = SymmetricTensor{2,3,Float64}((ε_2d[1, 1], ε_2d[1, 2], 0.0,
                        ε_2d[2, 2], 0.0, 0.0));
                else
                    ε̄u = shape_symmetric_gradient(cellvalues_u, q_point, j);
                end
                Ke[i, j] += δε ⊡ Dᵉ ⊡ ε̄u * dΩᵤ;
            end
            Re[i] += δε ⊡ σ * dΩᵤ;
        end
    end
    symmetrize_lower!(Ke)
    return Ke, Re
end

function assemble_residual_u!(Re::Vector, cellvalues_u,ue::Vector, material::Material, state, state_old, store::Bool)
    nbase_u = getnbasefunctions(cellvalues_u);
    for q_point in 1:getnquadpoints(cellvalues_u)
        dΩᵤ = getdetJdV(cellvalues_u, q_point);
        ϕ = state[q_point].ϕ;
        if material.dim == 2
            ε_PlaneStrain = function_symmetric_gradient(cellvalues_u, q_point, ue);
            ∇ˢu = SymmetricTensor{2,3,Float64}((ε_PlaneStrain[1, 1], ε_PlaneStrain[1, 2], 0.0, ε_PlaneStrain[2, 2], 0.0, 0.0));
        elseif material.dim == 3
            ∇ˢu = function_symmetric_gradient(cellvalues_u, q_point, ue);
        else
            error("Invalid element dimension");
        end
        
        gd, _, _ = Degradation(ϕ, material);
        Ψ⁺, _, σ, _ = Constitutive(∇ˢu, material, gd);
        H = max(Ψ⁺, state_old[q_point].H);
        if store
            state[q_point] = HistoryVariable(H, ϕ);
        end
        for i in 1:nbase_u
            if material.dim == 2
                δε_2d = shape_symmetric_gradient(cellvalues_u, q_point, i);
                δε = SymmetricTensor{2,3,Float64}((δε_2d[1, 1], δε_2d[1, 2], 0.0, δε_2d[2, 2], 0.0, 0.0));
            else
                δε = shape_symmetric_gradient(cellvalues_u, q_point, i);
            end
            Re[i] += (δε ⊡ σ) * dΩᵤ;
        end
    end
    return Re
end

function assemble_element_ϕ!(Ke::Matrix, Re::Vector, cellvalues_ϕ, ϕe::Vector, material::Material, state, state_old)
    nbase_ϕ = getnbasefunctions(cellvalues_ϕ);
    Gc = material.Gc;
    ℓ = material.ℓ;
    for q_point in 1:getnquadpoints(cellvalues_ϕ)
        dΩᵩ = getdetJdV(cellvalues_ϕ, q_point);
        ϕ = function_value(cellvalues_ϕ, q_point, ϕe);
        ∇ϕ = function_gradient(cellvalues_ϕ, q_point, ϕe);

        H = state[q_point].H;
        state[q_point] = HistoryVariable(H, ϕ);
        _, gd′, gd′′ = Degradation(ϕ, material);
        for i in 1:nbase_ϕ
            δϕ = shape_value(cellvalues_ϕ, q_point, i);
            ∇δϕ = shape_gradient(cellvalues_ϕ, q_point, i);
            for j in 1:i
                ϕ′ = shape_value(cellvalues_ϕ, q_point, j);
                ∇ϕ′ = shape_gradient(cellvalues_ϕ, q_point, j);
                Ke[i, j] += (gd′′ * ϕ′ * H * δϕ + Gc / ℓ * δϕ * ϕ′ + Gc * ℓ * ∇δϕ ⋅ ∇ϕ′) * dΩᵩ;
            end
            Re[i] += (gd′ * H * δϕ + Gc / ℓ * δϕ * ϕ + Gc * ℓ * ∇δϕ ⋅ ∇ϕ) * dΩᵩ;
        end
    end
    symmetrize_lower!(Ke)
    return Ke, Re
end

function assemble_residual_ϕ!(Re::Vector, cellvalues_ϕ, ϕe::Vector, material::Material, state, state_old, store::Bool)
    nbase_ϕ = getnbasefunctions(cellvalues_ϕ);
    Gc = material.Gc;
    ℓ = material.ℓ;
    for q_point in 1:getnquadpoints(cellvalues_ϕ);
        dΩᵩ = getdetJdV(cellvalues_ϕ, q_point);
        ϕ = function_value(cellvalues_ϕ, q_point, ϕe);
        ∇ϕ = function_gradient(cellvalues_ϕ, q_point, ϕe);
        H = state[q_point].H;
        if store
            state[q_point] = HistoryVariable(H, ϕ);
        end       
        _, gd′, _ = Degradation(ϕ, material);
        for i in 1:nbase_ϕ
            δϕ = shape_value(cellvalues_ϕ, q_point, i);
            ∇δϕ = shape_gradient(cellvalues_ϕ, q_point, i);
            Re[i] += (gd′ * H * δϕ +  Gc / ℓ * δϕ * ϕ +  Gc * ℓ * ∇δϕ ⋅ ∇ϕ) * dΩᵩ;
        end
    end
    return Re
end

function symmetrize_lower!(K)
    for i in 1:size(K,1)
        for j in i+1:size(K,1)
            K[i,j] = K[j,i]
        end
    end
end


function assemble_global_r(q::Vector, cellvalues, dh::DofHandler, material::Material, states, states_old, tag::String, store::Bool=true)
    nbase = getnbasefunctions(cellvalues);
    Re = zeros(nbase);
    R = zeros(ndofs(dh));
    for (i, cell) in enumerate(CellIterator(dh))
        reinit!(cellvalues, cell);
        fill!(Re, 0);
        eldofs = celldofs(cell);
        qe = q[eldofs];
        state = @view states[:, i];
        state_old = @view states_old[:, i];
        if tag == "u"
            assemble_residual_u!(Re, cellvalues, qe, material, state, state_old, store);
        elseif tag == "ϕ"
            assemble_residual_ϕ!(Re, cellvalues, qe, material, state, state_old, store);
        end
        R[eldofs] += Re;
    end
    return R
end;

function CreatCellvalues(ElementShape, ElementOrder, QuadratureOrder, dim)
	ip = Lagrange{dim, ElementShape, ElementOrder}()
	qr = QuadratureRule{dim, ElementShape}(QuadratureOrder)
	cellvalues_u = CellVectorValues(qr, ip)
	cellvalues_ϕ = CellScalarValues(qr, ip)
	return ip, qr, cellvalues_u, cellvalues_ϕ
end

function CreatBC(grid, dim)
	dbc₁ = Dirichlet(:u, getfaceset(grid, "top"), (x, t) -> t, 2)
	dbc₂ = Dirichlet(:u, getfaceset(grid, "bottom"), (x, t) -> [0, 0], [1, 2])

	dh_u = DofHandler(grid)
	add!(dh_u, :u, dim) 
	close!(dh_u)

	bch_u = ConstraintHandler(dh_u)
	add!(bch_u, dbc₁)
	add!(bch_u, dbc₂) 
	close!(bch_u)
	update!(bch_u, 0.0)

	dh_ϕ = DofHandler(grid)
	add!(dh_ϕ, :ϕ, 1) 
	close!(dh_ϕ)

	bch_ϕ = ConstraintHandler(dh_ϕ) 
	close!(bch_ϕ)
	update!(bch_ϕ, 0.0)

    return dh_u, dh_ϕ, bch_u, bch_ϕ
end

function Degradation(ϕ, mat)
    # kmin = 1e-15;
    kmin = 1e-8;
    flagD = mat.flagD
    if flagD == QuadraticDegradation
        gd = (1.0 - ϕ)^2 + kmin; 
        gd′ = -2.0(1.0 - ϕ);
        gd′′ = 2.0;
    elseif flagD == WuDegradation
        a₁ = mat.a₁;
        fact = (1.0 - ϕ)^2 + a₁*ϕ*(1 - 0.5ϕ);
        gd = (1.0 - ϕ)^2/fact + kmin;
        gd′ = -a₁*(1.0 - ϕ)/fact^2;
        gd′′ = 2a₁^2/fact^3 - 3a₁/fact^2;
        gd′′ = gd′′>0 ? gd′′ : 0;
    end
    return gd, gd′, gd′′
end

function Constitutive(ε::SymmetricTensor{2,3,Float64}, mat, gdn)
    Heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x, 0.5)))
    flag = mat.flag
    if flag == Isotropic
        μ = mat.G;
        λ = mat.λ;
        I = one(SymmetricTensor{2,3});
        𝕀 = one(SymmetricTensor{4,3});
        D = λ * I ⊗ I + 2μ * 𝕀;
        Ψ⁺ = 0.5 * ε ⊡ D ⊡ ε;
        Dᵉ = gdn * D; 
        σ = Dᵉ ⊡ ε;
        σ⁺ = σ;
    elseif flag == VolDev 
        K = mat.K;
        G = mat.G;
        I = one(SymmetricTensor{2,3});
        𝕀 = one(SymmetricTensor{4,3});
        D⁺ = K * Heaviside(tr(ε)) * I ⊗ I + 2G * (𝕀 - 1 / 3 * I ⊗ I);
        D⁻ = K * Heaviside(-tr(ε)) * I ⊗ I;
        Ψ⁺ = tr(ε) >0 ? 0.5*K*tr(ε)^2 + G*dev(ε) ⊡ dev(ε) : G*dev(ε) ⊡ dev(ε);
        σ⁺ = tr(ε) >0 ? K*tr(ε)*I + 2G*dev(ε) : 2G*dev(ε); 
        σ⁻ = tr(ε) <0 ? K*tr(ε)*I  : zero(Tensor{2, 3}); 
        Dᵉ = gdn * D⁺ + D⁻;
        σ = gdn *σ⁺ + σ⁻;
    elseif flag == Spectral 
        εₙ, Vₙ = eigen(ε);        # εₙ = eigvals(ε)
        μ = mat.G;
        λ = mat.λ;
        bracket₊(a::AbstractFloat) = a > 0 ? a : 0;
        bracket₋(a::AbstractFloat) = a < 0 ? a : 0;
        H₁₂⁺(x::AbstractFloat, y::AbstractFloat) = x ≠ y ? (bracket₊(x) - bracket₊(y)) / (x - y) : Heaviside(x);
        H₁₂⁻(x::AbstractFloat, y::AbstractFloat) = x ≠ y ? (bracket₋(x) - bracket₋(y)) / (x - y) : Heaviside(-x);
        I = one(SymmetricTensor{2,3});
        Ψ⁺ = tr(ε) >0 ? λ/2*tr(ε)^2 : 0.0;
        for e in εₙ;Ψ⁺ += e>0 ? μ*e^2 : 0.0;end;
        σ⁺ = λ*bracket₊(tr(ε))*I + 2μ*(bracket₊(εₙ[1])*Vₙ[:,1]⊗Vₙ[:,1] + bracket₊(εₙ[2])*Vₙ[:,2]⊗Vₙ[:,2] + bracket₊(εₙ[3])*Vₙ[:,3]⊗Vₙ[:,3]);
        σ⁻ = λ*bracket₋(tr(ε))*I + 2μ*(bracket₋(εₙ[1])*Vₙ[:,1]⊗Vₙ[:,1] + bracket₋(εₙ[2])*Vₙ[:,2]⊗Vₙ[:,2] + bracket₋(εₙ[3])*Vₙ[:,3]⊗Vₙ[:,3]);
        D⁺ = λ * Heaviside(tr(ε)) * I ⊗ I;
        D⁻ = λ * Heaviside(-tr(ε)) * I ⊗ I;
        for a in 1:3
            D⁺ += 2μ * Heaviside(εₙ[a]) * (Vₙ[:, a] ⊗ Vₙ[:, a]) ⊗ (Vₙ[:, a] ⊗ Vₙ[:, a]);
            D⁻ += 2μ * Heaviside(-εₙ[a]) * (Vₙ[:, a] ⊗ Vₙ[:, a]) ⊗ (Vₙ[:, a] ⊗ Vₙ[:, a]);
            for b in 1:3
                if a ≠ b
                    D⁺ += μ * H₁₂⁺(εₙ[a], εₙ[b]) * ((Vₙ[:, a] ⊗ Vₙ[:, b]) ⊗ (Vₙ[:, a] ⊗ Vₙ[:, b]) + (Vₙ[:, a] ⊗ Vₙ[:, b]) ⊗ (Vₙ[:, b] ⊗ Vₙ[:, a]));
                    D⁻ += μ * H₁₂⁻(εₙ[a], εₙ[b]) * ((Vₙ[:, a] ⊗ Vₙ[:, b]) ⊗ (Vₙ[:, a] ⊗ Vₙ[:, b]) + (Vₙ[:, a] ⊗ Vₙ[:, b]) ⊗ (Vₙ[:, b] ⊗ Vₙ[:, a]));
                end
            end
        end
        Dᵉ = gdn * D⁺ + D⁻;
        σ = gdn *σ⁺ + σ⁻;
    end
    return Ψ⁺, Dᵉ, σ, σ⁺
end

function Newton_raphson!(q::Vector, K::SparseMatrixCSC, cellvalues, dh::DofHandler, ch::ConstraintHandler, grid::Grid, Mat::Material, states, states_old, Ru_first, Rϕ_first, tag::String)
    iterations = 0;
    for nitr = 1:(solver.nitr_inner+1);
        if nitr > solver.nitr_inner;
            error("Reached maximum Newton iterations, aborting");
            break;
        end;
        K, r = assemble_global(q, cellvalues, K, dh, Mat, states, states_old, tag::String);
        apply_zero!(K, r, ch);
        norm_r = norm(r[Ferrite.free_dofs(ch)]);  # norm_r = maximum(abs.(r[Ferrite.free_dofs(ch)]))
        if tag =="u"
            TOL = solver.TOL_u * Ru_first;
        elseif tag =="ϕ"
            TOL = solver.TOL_ϕ * Rϕ_first;
        end
        if (norm_r < TOL) && (nitr > 1);
            break;
        end
        iterations += 1;
        K_active = K[Ferrite.free_dofs(ch),Ferrite.free_dofs(ch)];
		r_active =  r[Ferrite.free_dofs(ch)];
        # Δq = K_active \ r_active;
        q[Ferrite.free_dofs(ch)] -= K_active \ r_active;
    end
        print(tag*" converged in $iterations iterations \n");
    return q, iterations, K
end

function assemble_global(q::Vector, cellvalues, K::SparseMatrixCSC, dh::DofHandler, material::Material, states, states_old, tag::String)
    nbase = getnbasefunctions(cellvalues);

    Ke = zeros(nbase, nbase);
    Re = zeros(nbase);
    R = zeros(ndofs(dh));
    assembler = start_assemble(K, R);

    for (i, cell) in enumerate(CellIterator(dh))
        reinit!(cellvalues, cell);
        fill!(Ke, 0);
        fill!(Re, 0);
        eldofs = celldofs(cell);
        qe = q[eldofs];
        state = @view states[:, i];
        state_old = @view states_old[:, i];
        if tag == "u"
            assemble_element_u!(Ke, Re, cellvalues, qe, material, state, state_old);
        elseif  tag == "ϕ"
            assemble_element_ϕ!(Ke, Re, cellvalues, qe, material, state, state_old);
        end
        assemble!(assembler, eldofs, Ke, Re);
        
    end
    return K, R
end;

function Masterstep_u!(u, ϕ, K, R_u, R_ϕ, cellvalues_u, cellvalues_ϕ, dh, dh_u, dh_ϕ, ch, ch_u, ch_ϕ, material, states, states_old, flag_eta)
    @inbounds begin
        u_free = free_dofs(ch_u)
        ϕ_free = free_dofs(ch_ϕ)
        n_u = ndofs(dh_u)
        q = vcat(u, ϕ) 
        K = assemble_global_K_old(q, K, cellvalues_u, cellvalues_ϕ, dh, material, states, states_old)

        Kuu = @view K[1:n_u, 1:n_u]
        Kuϕ = @view K[1:n_u, n_u+1:end]
        Kϕu = @view K[n_u+1:end, 1:n_u]
        Kϕϕ = @view K[n_u+1:end, n_u+1:end]
        K_uu = @view (Kuu[u_free, u_free])
        K_uϕ = @view (Kuϕ[u_free, ϕ_free])
        K_ϕu = @view (Kϕu[ϕ_free, u_free])
        K_ϕϕ = @view (Kϕϕ[ϕ_free, ϕ_free])

        K_uu = sparse(Symmetric(K_uu))
        K_ϕϕ = sparse(Symmetric(K_ϕϕ))
        K_uu_fact = cholesky(K_uu)
        K_ϕϕ_fact = cholesky(K_ϕϕ)
        eta = 1;
        if flag_eta == 1
            rho = cal_rho1(K_uu, K_ϕϕ, K_uϕ, K_ϕu)
            if rho < 1
                println("eta = ", eta)
            elseif rho > 1
                eta = 0.9/sqrt(rho);
                println("eta = ", eta)
            end
        end
        y_active = -(K_uu_fact \ (R_u - eta * K_uϕ * (K_ϕϕ_fact \ R_ϕ)))
        Δu_active = gmres_user(K_uu_fact, K_ϕϕ_fact, eta * K_uϕ, eta * K_ϕu, y_active; max_iter=50, tol=1e-6)
        Δu = zeros(ndofs(dh_u))
        Δu[u_free] .= Δu_active   
    end
    return Δu
end

function assemble_coupled_K(u, ϕ, cellvalues_u, cellvalues_ϕ, dh_u, dh_ϕ, material, states, states_old)
    # 计算总自由度：位移 (u) + 相场 (ϕ)
    nu = ndofs(dh_u)
    nϕ = ndofs(dh_ϕ)
    n_dofs = nu + nϕ
    
    # 开始组装，不显式指定稀疏模式
    assembler = start_assemble()
    
    # 元素级的组装
    nbase_u = getnbasefunctions(cellvalues_u)  # 位移场的基函数数量
    nbase_ϕ = getnbasefunctions(cellvalues_ϕ)  # 相场基函数数量
    Ke = BlockArray(zeros(nbase_u + nbase_ϕ, nbase_u + nbase_ϕ), [nbase_u, nbase_ϕ], [nbase_u, nbase_ϕ])
    
    # 遍历单元
    for (i, cell) in enumerate(CellIterator(dh_u.grid))
        fill!(Ke, 0)  # 清空元素刚度矩阵
        reinit!(cellvalues_u, cell)  # 重置单元值
        reinit!(cellvalues_ϕ, cell)
        
        # 获取单元自由度
        eldofs_u = celldofs(dh_u, i)
        eldofs_ϕ = celldofs(dh_ϕ, i)
        eldofs = vcat(eldofs_u, eldofs_ϕ)  # 合并自由度
        
        # 获取单元解向量
        ue = u[eldofs_u]
        ϕe = ϕ[eldofs_ϕ]
        qe = vcat(ue, ϕe)  # 合并解向量
        
        state = @view states[:, i]
        state_old = @view states_old[:, i]
        
        # 组装元素刚度矩阵
        assemble_element_ke_old!(Ke, cell, cellvalues_u, cellvalues_ϕ, qe, material, state, state_old)
        
        # 将元素矩阵组装到全局矩阵
        assemble!(assembler, eldofs, Ke)
    end
    
    # 完成组装并返回全局矩阵
    K, _ = end_assemble(assembler)
    return K
end

function assemble_element_ke_old!(Ke, cell, cellvalues_u, cellvalues_ϕ, qe, material, state, state_old)
    nbase_u = getnbasefunctions(cellvalues_u)
    nbase_ϕ = getnbasefunctions(cellvalues_ϕ)
    u▄, ϕ▄ = 1, 2 #位移是块1，相场是块2
    reinit!(cellvalues_u, cell)
    reinit!(cellvalues_ϕ, cell)

    ϕe = qe[nbase_u+1:end]
    ue = qe[1:nbase_u]

    Gc = material.Gc
    ℓ = material.ℓ

    for q_point in 1:getnquadpoints(cellvalues_u)
        dΩᵤ = getdetJdV(cellvalues_u, q_point)
        if material.dim == 2
            ε_PlaneStrain = function_symmetric_gradient(cellvalues_u, q_point, ue)
            ∇ˢu = SymmetricTensor{2,3,Float64}((ε_PlaneStrain[1, 1], ε_PlaneStrain[1, 2], 0.0,ε_PlaneStrain[2, 2], 0.0, 0.0))
        elseif material.dim == 3
            ∇ˢu = function_symmetric_gradient(cellvalues_u, q_point, ue)
        else
            error("Invalid element dimension")
        end

        ϕ = state[q_point].ϕ
        gd, gd′, gd′′  = Degradation(ϕ, material);
        Ψ⁺, Dᵉ, σ, σ⁺ = Constitutive(∇ˢu, material, gd);
        H = max(Ψ⁺, state_old[q_point].H);
        #state[q_point] = HistoryVariable(H, ϕ);

        for i in 1:nbase_u
            if material.dim == 2
                δε_2d = shape_symmetric_gradient(cellvalues_u, q_point, i);
                δε = SymmetricTensor{2,3,Float64}((δε_2d[1, 1], δε_2d[1, 2], 0.0,δε_2d[2, 2], 0.0, 0.0));
            else
                δε = shape_symmetric_gradient(cellvalues_u, q_point, i);
            end
            for j in 1:nbase_u
                if material.dim == 2
                    ε_2d = shape_symmetric_gradient(cellvalues_u, q_point, j);
                    ε̄u = SymmetricTensor{2,3,Float64}((ε_2d[1, 1], ε_2d[1, 2], 0.0,
                        ε_2d[2, 2], 0.0, 0.0));
                else
                    ε̄u = shape_symmetric_gradient(cellvalues_u, q_point, j);
                end
                Ke[BlockIndex((u▄, u▄), (i, j))] += δε ⊡ Dᵉ ⊡ ε̄u * dΩᵤ  #ε 是 B and δε 是 Bᵀ 故 BᵀDB
            end
            for j in 1:nbase_ϕ
                ϕ′ = shape_value(cellvalues_ϕ, q_point, j)
                Ke[BlockIndex((u▄, ϕ▄), (i, j))] += (gd′ * δε ⊡ σ⁺ * ϕ′) * dΩᵤ
            end
        end
        for i in 1:nbase_ϕ
            δϕ = shape_value(cellvalues_ϕ, q_point, i)
            ∇δϕ = shape_gradient(cellvalues_ϕ, q_point, i)
            for j in 1:nbase_ϕ
                ϕ′ = shape_value(cellvalues_ϕ, q_point, j)
                ∇ϕ′ = shape_gradient(cellvalues_ϕ, q_point, j)
                Ke[BlockIndex((ϕ▄, ϕ▄), (i, j))] += (gd′′ * ϕ′ * H * δϕ + Gc / ℓ *δϕ * ϕ′ +  Gc *ℓ * ∇δϕ ⋅ ∇ϕ′) * dΩᵤ
            end
            for j in 1:nbase_u
                if material.dim == 2
                    ε_2d = shape_symmetric_gradient(cellvalues_u, q_point, j);
                    ε̄u = SymmetricTensor{2,3,Float64}((ε_2d[1, 1], ε_2d[1, 2], 0.0,
                        ε_2d[2, 2], 0.0, 0.0));
                else
                    ε̄u = shape_symmetric_gradient(cellvalues_u, q_point, j);
                end
                if state[q_point].H > state_old[q_point].H
                    Ke[BlockIndex((ϕ▄, u▄), (i, j))] += (gd′ *δϕ * σ⁺ ⊡ ε̄u) * dΩᵤ
                else
                    Ke[BlockIndex((ϕ▄, u▄), (i, j))] += 0.0
                end   
            end
        end
    end    
end

function closure_M(K_uu, K_ϕϕ, K_uϕ, K_ϕu, x)
    result = x - K_uu\(K_uϕ*(K_ϕϕ\(K_ϕu*x)));
    return result
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

# 截断函数：保留前 k 个最大的元素
function truncate_vector(x, k)
    # 确保 k 是整数
    k = Int(round(k))
    
    # 获取绝对值最大的 k 个元素的索引
    idx = sortperm(abs.(x), rev=true)[1:k]
    
    # 创建一个新的向量，只保留前 k 个最大的元素，其余置零
    truncated_x = zeros(length(x))
    truncated_x[idx] = x[idx]
    return truncated_x
end


function gmres_user( K11, K22, K12, K21, b::Vector{Float64}; max_iter::Int=50, tol::Float64=1e-6)::Vector{Float64}
    n = length(b)
    x0 = zeros(Float64, n) 
    x_solution = copy(x0) 
    r0_vec = copy(b) 
    β = norm(r0_vec)

    V_basis = zeros(Float64, n, max_iter + 1)
    H_hessenberg = zeros(Float64, max_iter + 1, max_iter)
    g_rhs = zeros(Float64, max_iter + 1)    
    
    cs_givens = zeros(Float64, max_iter) 
    sn_givens = zeros(Float64, max_iter) 
    
    axpy!(1.0 / β, r0_vec, @view V_basis[:, 1]) # V_basis[:,1] = r0_vec / β
    g_rhs[1] = β

    j_actual_iters = 0 
    @inbounds begin
    for j in 1:max_iter 
        j_actual_iters = j
        current_V_j_view = @view V_basis[:, j]
        w = closure_M(K11, K22, K12, K21, current_V_j_view)

        # --- 正交化: Modified Gram-Schmidt (MGS) ---
        for i in 1:j
            V_i_view = @view V_basis[:,i]
            h_coeff = dot(V_i_view, w) 
            H_hessenberg[i, j] = h_coeff 
            axpy!(-h_coeff, V_i_view, w) 
        end
        # --- MGS 结束 ---
        
        H_hessenberg[j+1, j] = norm(w)

        if H_hessenberg[j+1, j] < eps(Float64) * β 
            break 
        end
        
        rmul!(w, 1.0 / H_hessenberg[j+1, j]) 
        (@view V_basis[:, j+1]) .= w   

        for i_rot in 1:(j-1)
            h_i_j_curr = H_hessenberg[i_rot, j]
            h_ip1_j_curr = H_hessenberg[i_rot+1, j]
            cs_val = cs_givens[i_rot]
            sn_val = sn_givens[i_rot]
            H_hessenberg[i_rot, j]   = cs_val * h_i_j_curr + sn_val * h_ip1_j_curr
            H_hessenberg[i_rot+1, j] = -sn_val * h_i_j_curr + cs_val * h_ip1_j_curr
        end
        
        a_rot = H_hessenberg[j,j]
        b_rot = H_hessenberg[j+1,j]
        
        local current_c, current_s
        rho_rot = sqrt(a_rot^2 + b_rot^2) 
        if rho_rot != 0.0 
            current_c = a_rot / rho_rot
            current_s = b_rot / rho_rot
            H_hessenberg[j,j] = rho_rot
            H_hessenberg[j+1,j] = 0.0 
        else
            current_c = 1.0
            current_s = 0.0
        end
        cs_givens[j] = current_c
        sn_givens[j] = current_s
        
        g_rhs[j+1] = -current_s * g_rhs[j] 
        g_rhs[j]   =  current_c * g_rhs[j] 

        relres = abs(g_rhs[j+1]) / β

        if relres < tol
            println("第 $j 次迭代收敛，相对残差为 $relres")
            break
        end
    end
    end # @inbounds end
    if j_actual_iters == 0
        if β < tol return x_solution end
        return x_solution 
    end
    
    H_solve_view = @view H_hessenberg[1:j_actual_iters, 1:j_actual_iters]
    g_solve_view = @view g_rhs[1:j_actual_iters]
    
    y_solution = UpperTriangular(H_solve_view) \ g_solve_view

    mul!(x_solution, (@view V_basis[:, 1:j_actual_iters]), y_solution, 1.0, 0.0)
    return x_solution
end

function assemble_global_K_old(q, K, cellvalues_u,cellvalues_ϕ, dh, material, states, states_old)

    assembler = start_assemble(K)
    nu = getnbasefunctions(cellvalues_u)
    nϕ = getnbasefunctions(cellvalues_ϕ)

    Ke = BlockArray(zeros(nu + nϕ, nu + nϕ), [nu, nϕ], [nu, nϕ]) # local stiffness matrix
    for (i, cell)  in enumerate(CellIterator(dh))
        fill!(Ke, 0)
        reinit!(cellvalues_ϕ, cell);
        reinit!(cellvalues_u, cell);
        eldofs = celldofs(cell)
        qe = q[eldofs]
        state = @view states[:, i]
        state_old = @view states_old[:, i]
        assemble_element_ke_old!(Ke, cell, cellvalues_u, cellvalues_ϕ, qe, material, state, state_old)

        assemble!(assembler, eldofs, Ke)
    end
    return K
end

function OutputForce(u, cellvalues_u, dh_u, grid, material, states, states_old, set::String, tag)
	F = assemble_global_r(u, cellvalues_u, dh_u, material, states, states_old, tag, true);
	F_x = 0.0;
	F_y = 0.0;
	Fnodal = reshape_to_nodes(dh_u, F, :u);
    # Fnodal = evaluate_at_grid_nodes(dh, F, :u);
	if set ∈ keys(grid.nodesets);
		outputset = grid.nodesets[set];
	elseif set ∈ keys(grid.facesets[set]);
		print("facesets are currently not supported for force output");
	else;
		print("Warning invalid set for force output");
	end;
	for (i, n) in enumerate(grid.nodesets[set]);
		F_x += Fnodal[1, n];
		F_y += Fnodal[2, n];
		# F_x += Fnodal[i][1];
		# F_y += Fnodal[i][2];
	end;
	return F_x, F_y;
end;
function CrackTrack(q, dh, cellValues, grid, a0, CrackDir)
	Ac = a0
	v = CrackDir
	for (i, cell) in enumerate(CellIterator(dh))
		reinit!(cellValues, cell)
		eldofs = celldofs(cell)
		ϕe = q[eldofs]
		if maximum(ϕe) >= 0.95
			node_coords = getcoordinates(grid, i)
			for q_point in 1:getnquadpoints(cellValues)
				#Phase field value
				ϕ = function_value(cellValues, q_point, ϕe)
				if ϕ >= 0.95
					coords = spatial_coordinate(cellValues, q_point, node_coords)
					Ac = coords[v] > Ac ? coords[v] : Ac
				end
			end
		end
	end
	return Ac
end
function write_to_txt(filename, variable)
    open(filename, "w") do file
        write(file, string(variable))  # 将变量转换为字符串写入文件
    end
end
function save_results(folder, timestep, u, ϕ, states, states_old)
    # 确保文件夹存在，如果不存在则创建
    if !isdir(folder)
        mkdir(folder)
    end

    # 将时间步格式化为字符串，保留 5 位小数，并用下划线替代小数点
    timestep_str = replace(@sprintf("%.5f", timestep), "." => "_")
    filename = joinpath(folder, "time_$timestep_str.jld")  # 将文件保存到指定文件夹
    save(filename, "timestep", timestep, "u", u, "ϕ", ϕ, "states", states, "states_old", states_old)
end
function load_results(folder, timestep)
    # 将时间步格式化为字符串，保留 5 位小数，并用下划线替代小数点
    timestep_str = replace(@sprintf("%.5f", timestep), "." => "_")
    filename = joinpath(folder, "time_$timestep_str.jld")  # 从指定文件夹加载文件
    if !isfile(filename)
        error("File $filename does not exist.")
    end
    data = load(filename)
    return data["timestep"], data["u"], data["ϕ"], data["states"], data["states_old"]
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
        if material.dim == 2
            ε_PlaneStrain = function_symmetric_gradient(cellvalues_u, q_point, ue);
            ∇ˢu = SymmetricTensor{2,3,Float64}((ε_PlaneStrain[1, 1], ε_PlaneStrain[1, 2], 0.0, ε_PlaneStrain[2, 2], 0.0, 0.0));
        elseif material.dim == 3
            ∇ˢu = function_symmetric_gradient(cellvalues_u, q_point, ue) - εθ*one(SymmetricTensor{2, 3});
        else
            error("Invalid element dimension");
        end
        
        gd, _, _ = Degradation(ϕ, material);
        Ψ⁺, _, σ, _ = Constitutive(∇ˢu, material, gd);
        H = max(Ψ⁺, state_old[q_point].H);
        state[q_point] = HistoryVariable(H, ϕ);
    end
    return state
end
