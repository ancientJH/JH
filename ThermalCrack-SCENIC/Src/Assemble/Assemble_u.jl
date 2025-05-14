function assemble_element_u!(Ke::Matrix, Re::Vector, cellvalues_u,ue::Vector, material::Material, state, state_old)
    nbase_u = getnbasefunctions(cellvalues_u);
    for q_point in 1:getnquadpoints(cellvalues_u)
        dΩᵤ = getdetJdV(cellvalues_u, q_point);
        ϕ = state[q_point].ϕ;
        θ = state[q_point].θ;
        εθ = material.α*(θ - material.θ₀);
        if material.dim == 2
            ε_PlaneStrain = function_symmetric_gradient(cellvalues_u, q_point, ue);
            ε_PlaneStrain -= εθ*one(SymmetricTensor{2, 2});
            ∇ˢu = SymmetricTensor{2,3,Float64}((ε_PlaneStrain[1, 1], ε_PlaneStrain[1, 2], 0.0, ε_PlaneStrain[2, 2], 0.0, 0.0));
        elseif material.dim == 3
            ∇ˢu = function_symmetric_gradient(cellvalues_u, q_point, ue);
            ∇ˢu -= εθ*one(SymmetricTensor{2, 3});
        else
            error("Invalid element dimension");
        end
        gd, _, _ = Degradation(ϕ, material);
        Ψ⁺, Dᵉ, σ ,_ = Constitutive(∇ˢu, material, gd);
        H = max(Ψ⁺, state_old[q_point].H);
        state[q_point] = HistoryVariable(H, ϕ, θ);

        for i in 1:nbase_u
            if material.dim == 2
                δε_2d = shape_symmetric_gradient(cellvalues_u, q_point, i);
                δεu = SymmetricTensor{2,3,Float64}((δε_2d[1, 1], δε_2d[1, 2], 0.0,
                δε_2d[2, 2], 0.0, 0.0));
            else
                δεu = shape_symmetric_gradient(cellvalues_u, q_point, i);
            end
            for j in 1:nbase_u
                if material.dim == 2
                    ε_2d = shape_symmetric_gradient(cellvalues_u, q_point, j);
                    ε̄u = SymmetricTensor{2,3,Float64}((ε_2d[1, 1], ε_2d[1, 2], 0.0,
                        ε_2d[2, 2], 0.0, 0.0));
                else
                    ε̄u = shape_symmetric_gradient(cellvalues_u, q_point, j);
                end
                Ke[i, j] += δεu ⊡ Dᵉ ⊡ ε̄u * dΩᵤ;
            end
            Re[i] += δεu ⊡ σ * dΩᵤ;
        end
    end
    # symmetrize_lower!(Ke)
    return Ke, Re
end

function assemble_residual_u!(Re::Vector, cellvalues_u,ue::Vector, material::Material, state, state_old, store::Bool)
    nbase_u = getnbasefunctions(cellvalues_u);
    for q_point in 1:getnquadpoints(cellvalues_u)
        dΩᵤ = getdetJdV(cellvalues_u, q_point);
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
        if store
            state[q_point] = HistoryVariable(H, ϕ, θ);
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
