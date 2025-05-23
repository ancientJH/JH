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

function assemble_global_Ku(u, K, cellvalues_u, dh_u, material, states, states_old)

    assembler = start_assemble(K)
    nu = getnbasefunctions(cellvalues_u)

    Ke = zeros(nu, nu) # local stiffness matrix
    for (i, cell)  in enumerate(CellIterator(dh_u))
        fill!(Ke, 0)
        reinit!(cellvalues_u, cell);
        eldofs = celldofs(dh_u, i)
        ue = u[eldofs]
        state = @view states[:, i]
        state_old = @view states_old[:, i]
        assemble_element_Ku!(Ke, cellvalues_u, ue, material, state, state_old)
        assemble!(assembler, eldofs, Ke)
    end
    return K
end

function assemble_global_Kϕ(ϕ, K, cellvalues_ϕ, dh_ϕ, material, states, states_old)
    assembler = start_assemble(K)
    nϕ = getnbasefunctions(cellvalues_ϕ)
    Ke = zeros(nϕ, nϕ) # local stiffness matrix
    for (i, cell)  in enumerate(CellIterator(dh_ϕ))
        fill!(Ke, 0)
        reinit!(cellvalues_ϕ, cell);
        eldofs = celldofs(dh_ϕ, i)
        ϕe = ϕ[eldofs]
        state = @view states[:, i]
        state_old = @view states_old[:, i]
        assemble_element_Kϕ!(Ke, cellvalues_ϕ, ϕe, material, state, state_old)
        assemble!(assembler, eldofs, Ke)
    end
    return K
end