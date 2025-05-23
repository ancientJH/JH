
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

function assemble_residual_ϕ_withu!(Re::Vector, cellvalues_u, cellvalues_ϕ, ue::Vector, ϕe::Vector, material::Material, state, state_old, store::Bool)
    nbase_ϕ = getnbasefunctions(cellvalues_ϕ);
    Gc = material.Gc;
    ℓ = material.ℓ;
    @inbounds for q_point in 1:getnquadpoints(cellvalues_ϕ);
        dΩᵩ = getdetJdV(cellvalues_ϕ, q_point);
        ϕ = function_value(cellvalues_ϕ, q_point, ϕe);
        ∇ϕ = function_gradient(cellvalues_ϕ, q_point, ϕe);
        gd, gd′, _ = Degradation(ϕ, material);
        if material.dim == 2
            ε_PlaneStrain = function_symmetric_gradient(cellvalues_u, q_point, ue) 
            ∇ˢu = SymmetricTensor{2,3,Float64}((ε_PlaneStrain[1, 1], ε_PlaneStrain[1, 2], 0.0, ε_PlaneStrain[2, 2], 0.0, 0.0));
        elseif material.dim == 3
            ∇ˢu = function_symmetric_gradient(cellvalues_u, q_point, ue);
        else
            error("Invalid element dimension");
        end
        Ψ⁺, _, _ = Constitutive(∇ˢu, material, gd);
        H = max(Ψ⁺, state_old[q_point].H);
        if store
            state[q_point] = HistoryVariable(H, ϕ);
        end       
        for i in 1:nbase_ϕ
            δϕ = shape_value(cellvalues_ϕ, q_point, i);
            ∇δϕ = shape_gradient(cellvalues_ϕ, q_point, i);
            Re[i] += (gd′ * H * δϕ +  Gc / ℓ * δϕ * ϕ +  Gc * ℓ * ∇δϕ ⋅ ∇ϕ) * dΩᵩ;
        end
    end
    return Re
end

function assemble_element_Kϕ!(Ke::Matrix, cellvalues_ϕ, ϕe::Vector, material::Material, state, state_old)
    nbase_ϕ = getnbasefunctions(cellvalues_ϕ);
    Gc = material.Gc;
    ℓ = material.ℓ;
    @inbounds for q_point in 1:getnquadpoints(cellvalues_ϕ)
        dΩᵩ = getdetJdV(cellvalues_ϕ, q_point);
        ϕ = function_value(cellvalues_ϕ, q_point, ϕe);

        H = state[q_point].H;
        state[q_point] = HistoryVariable(H, ϕ);
        _, _, gd′′ = Degradation(ϕ, material);
        for i in 1:nbase_ϕ
            δϕ = shape_value(cellvalues_ϕ, q_point, i);
            ∇δϕ = shape_gradient(cellvalues_ϕ, q_point, i);
            for j in 1:i
                ϕ′ = shape_value(cellvalues_ϕ, q_point, j);
                ∇ϕ′ = shape_gradient(cellvalues_ϕ, q_point, j);
                Ke[i, j] += (gd′′ * ϕ′ * H * δϕ + Gc / ℓ * δϕ * ϕ′ + Gc * ℓ * ∇δϕ ⋅ ∇ϕ′) * dΩᵩ;
            end
        end
    end
    symmetrize_lower!(Ke)
    return Ke
end