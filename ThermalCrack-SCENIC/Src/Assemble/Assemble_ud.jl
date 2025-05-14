function assemble_element_ke_ud!(Ke, cell, cellvalues_u, cellvalues_ϕ, cellvalues_θ, qe, material, state, state_old)
    nbase_u = getnbasefunctions(cellvalues_u)
    nbase_ϕ = getnbasefunctions(cellvalues_ϕ)
    #nbase_θ = getnbasefunctions(cellvalues_θ)
    u▄, ϕ▄ = 1, 2
    reinit!(cellvalues_u, cell)
    reinit!(cellvalues_ϕ, cell)
    reinit!(cellvalues_θ, cell)
    ϕe = qe[nbase_u + 1 : end]
    #θe = qe[nbase_u + nbase_ϕ + 1 : end]
    ue = qe[1 : nbase_u] 

    Gc = material.Gc
    ℓ = material.ℓ

    for q_point in 1:getnquadpoints(cellvalues_u)
        
        ϕ = state[q_point].ϕ;
        θ = state[q_point].θ;
        # ϕ = function_value(cellvalues_ϕ, q_point, ϕe);
        # θ = function_value(cellvalues_θ, q_point, θe);
        # ∇θ = function_gradient(cellvalues_θ, q_point, θe);

        εθ = material.α*(θ - material.θ₀);

        dΩᵤ = getdetJdV(cellvalues_u, q_point)
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


        gd, gd′, gd′′ = Degradation(ϕ, material);
        Ψ⁺, Dᵉ, σ, σ⁺ = Constitutive(∇ˢu, material, gd)
        H = max(Ψ⁺, state_old[q_point].H)

        state[q_point] = HistoryVariable(H, ϕ, θ);

        for i in 1:nbase_u
            if material.dim == 2
                δε_2d = shape_symmetric_gradient(cellvalues_u, q_point, i)
                δεu = SymmetricTensor{2,3,Float64}((δε_2d[1, 1], δε_2d[1, 2], 0.0,
                    δε_2d[2, 2], 0.0, 0.0))
            else
                δεu = shape_symmetric_gradient(cellvalues_u, q_point, i)
            end
            for j in 1:nbase_u
                if material.dim == 2
                    ε_2d = shape_symmetric_gradient(cellvalues_u, q_point, j)
                    ε̄u = SymmetricTensor{2,3,Float64}((ε_2d[1, 1], ε_2d[1, 2], 0.0,
                        ε_2d[2, 2], 0.0, 0.0))
                else
                    ε̄u = shape_symmetric_gradient(cellvalues_u, q_point, j)
                end
                Ke[BlockIndex((u▄, u▄), (i, j))] += δεu ⊡ Dᵉ ⊡ ε̄u * dΩᵤ  #ε 是 B and δε 是 Bᵀ 故 BᵀDB
            end
            for j in 1:nbase_ϕ
                ϕ′ = shape_value(cellvalues_ϕ, q_point, j)
                Ke[BlockIndex((u▄, ϕ▄), (i, j))] += (gd′ * δεu ⊡ σ⁺ * ϕ′) * dΩᵤ
            end
        end
        for i in 1:nbase_ϕ
            δϕ = shape_value(cellvalues_ϕ, q_point, i)
            ∇δϕ = shape_gradient(cellvalues_ϕ, q_point, i)
            for j in 1:nbase_ϕ
                ϕ′ = shape_value(cellvalues_ϕ, q_point, j)
                ∇ϕ′ = shape_gradient(cellvalues_ϕ, q_point, j)
                Ke[BlockIndex((ϕ▄, ϕ▄), (i, j))] += (gd′′ * ϕ′ * H * δϕ + Gc / ℓ *(δϕ * ϕ′ +  ℓ *ℓ * ∇δϕ ⋅ ∇ϕ′)) * dΩᵤ
            end
            for j in 1:nbase_u
                if material.dim == 2
                    ε_2d = shape_symmetric_gradient(cellvalues_u, q_point, j)
                    ε̄u = SymmetricTensor{2,3,Float64}((ε_2d[1, 1], ε_2d[1, 2], 0.0,
                        ε_2d[2, 2], 0.0, 0.0))
                else
                    ε̄u = shape_symmetric_gradient(cellvalues_u, q_point, j)
                end
                if state[q_point].H > state_old[q_point].H
                    Ke[BlockIndex((ϕ▄, u▄), (i, j))] += (gd′ *δϕ * σ⁺ ⊡ ε̄u) * dΩᵤ
                end   
            end
        end
    end    
end

function assemble_global_K_ud(q, K, cellvalues_u,cellvalues_ϕ, cellvalues_θ, dhud, material, states, states_old)
# 这里的q传的是u和d
    assembler = start_assemble(K)
    nu = getnbasefunctions(cellvalues_u)
    nϕ = getnbasefunctions(cellvalues_ϕ)
    # nθ = getnbasefunctions(cellvalues_θ)

    Ke = BlockArray(zeros(nu + nϕ , nu + nϕ), [nu, nϕ], [nu, nϕ]) # local stiffness matrix
    for (i, cell)  in enumerate(CellIterator(dhud))
        fill!(Ke, 0)
        eldofs = celldofs(cell)
        qe = q[eldofs]
        state = @view states[:, i]
        state_old = @view states_old[:, i]
        assemble_element_ke_ud!(Ke, cell, cellvalues_u, cellvalues_ϕ, cellvalues_θ, qe, material, state, state_old)

        assemble!(assembler, eldofs, Ke)
    end
    return K
end