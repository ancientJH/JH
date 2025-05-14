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
    ρ::T  # mass density
    Cp::T # specific heat capacity
    k₀::T # thermal conductivity
    α::T  # thermal expansion coefficient
    θ₀::T # initial tempetarture
    Δt::T # time step
    s::T  # parameter for Hughes degradation 
    a₁::T # parameter for Wu degradation 
    flag::StrainDecomp
    flagD::DegradType
    dim::Int64
    
end
function Material(E, ν, Gc, σc, ℓ, ρ, Cp, k₀, α, θ₀, Δt, s, flag, flagD, dim)
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)
    λ = K - 2G / 3
    a₁ = 27E*Gc/(128σc^2*ℓ₀);
    return Material(E, G, K, λ, ν, Gc, σc, ℓ, ρ, Cp, k₀, α, θ₀, Δt, s, a₁, flag, flagD, dim)
end

#定义历史变量结构体
struct HistoryVariable{T}
	H::T # History variable
	ϕ::T #phase field variable from last increment
	θ::T #temperature variable from last increment
end
function HistoryVariable()
	return HistoryVariable(0.0, 0.0, 0.0)
end

#定义求解参数结构体
mutable struct SolverState{T,F}
    loadsteps::Vector{F}
    nitr_inner::T
    nitr_outer::T
    TOL_u::F
    TOL_ϕ::F
    TOL_θ::F
    # Ndofϕ::Int64  #每个单元拥有的ϕ自由度
end
function SolverState(loadsteps, nitr_inner, nitr_outer, TOL_u, TOL_ϕ, TOL_θ)
    return SolverState(loadsteps, nitr_inner, nitr_outer, TOL_u, TOL_ϕ, TOL_θ)
end

#定义输出参数结构体
mutable struct OutputVariables{T}
    plotframe::T
    totalIterations_outer::T
    totalIterations_inter::T
    totalIterations_ϕ::T
    totalIterations_θ::T
    totalIterations_u::T
    plotFrequency::T
    historyFrequency::T
    a0::Float64
    CrackDir::T
    OutputSet::String
end
#输出参数结构体的初始化
function OutputVariables(field_frequency, history_frequency, a0, CrackDir, outputset)
    return OutputVariables(0, 0, 0, 0, 0, 0, field_frequency, history_frequency, a0, CrackDir, outputset)
end