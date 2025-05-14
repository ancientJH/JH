function assemble_element_θ!(Ke::Matrix, Re::Vector, cellvalues_θ, θe::Vector, material::Material, state, state_old)
    nbase_θ = getnbasefunctions(cellvalues_θ);
    ρ = material.ρ;
    Cp = material.Cp;
    Δt = material.Δt;
    k₀ = material.k₀;
    for q_point in 1:getnquadpoints(cellvalues_θ)
        dΩθ = getdetJdV(cellvalues_θ, q_point);
        θ = function_value(cellvalues_θ, q_point, θe);
        ∇θ = function_gradient(cellvalues_θ, q_point, θe);
        ϕ = state[q_point].ϕ;
        H = state[q_point].H;
        θₙ = state_old[q_point].θ;
        state[q_point] = HistoryVariable(H, ϕ, θ);
        gd, _, _ = Degradation(ϕ, material);
        for i in 1:nbase_θ
            δθ = shape_value(cellvalues_θ, q_point, i);
            ∇δθ = shape_gradient(cellvalues_θ, q_point, i);
            for j in 1:i
                θ′ = shape_value(cellvalues_θ, q_point, j);
                ∇θ′ = shape_gradient(cellvalues_θ, q_point, j);
                Ke[i, j] += (ρ * Cp/Δt * δθ * θ′  + gd * k₀ * ∇δθ ⋅ ∇θ′) * dΩθ;
            end
            Re[i] += (ρ * Cp * (θ - θₙ)/Δt * δθ +  gd * k₀* ∇δθ ⋅ ∇θ) * dΩθ;
        end
    end
    symmetrize_lower!(Ke)
    return Ke, Re
end

function assemble_residual_θ!(Re::Vector, cellvalues_θ, θe::Vector, material::Material, state, state_old, store::Bool)
    nbase_θ = getnbasefunctions(cellvalues_θ);
    ρ = material.ρ;
    Cp = material.Cp;
    Δt = material.Δt;
    k₀ = material.k₀;
    for q_point in 1:getnquadpoints(cellvalues_θ);
        dΩθ = getdetJdV(cellvalues_θ, q_point);
        θ = function_value(cellvalues_θ, q_point, θe);
        ∇θ = function_gradient(cellvalues_θ, q_point, θe);
        ϕ = state[q_point].ϕ;
        H = state[q_point].H;
        θₙ = state_old[q_point].θ;
        if store
            state[q_point] = HistoryVariable(H, ϕ, θ);
        end       
        gd, _, _ = Degradation(ϕ, material);
        for i in 1:nbase_θ
            δθ = shape_value(cellvalues_θ, q_point, i);
            ∇δθ = shape_gradient(cellvalues_θ, q_point, i);
            Re[i] += (ρ * Cp * (θ - θₙ)/Δt * δθ +  gd * k₀* ∇δθ ⋅ ∇θ) * dΩθ;
        end
    end
    return Re
end