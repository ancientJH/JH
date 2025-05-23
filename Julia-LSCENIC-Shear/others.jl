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

function updateϕ(cellvalues_ϕ, ϕ, dh_ϕ, material, states, states_old)
    for (i, cell) in enumerate(CellIterator(dh_ϕ))
        reinit!(cellvalues_ϕ, cell)
        eldofs_ϕ = celldofs(cell)         # 针对dh_u，编号就是u场的本地编号
        ϕe = ϕ[eldofs_ϕ]
        state = @view states[:, i]
        state_old = @view states_old[:, i]
        states[:, i] = assemble_element_hisϕ!(cellvalues_ϕ, ϕe, material, state, state_old)
    end
    return states
end

function assemble_element_hisϕ!(cellvalues_ϕ,ϕe::Vector, material::Material, state, state_old)
    for q_point in 1:getnquadpoints(cellvalues_ϕ)
        ϕ = function_value(cellvalues_ϕ, q_point, ϕe);
        H = state[q_point].H;
        state[q_point] = HistoryVariable(H, ϕ);
    end
    return state
end