function assemble_element_ke!(Ke, cell, cellvalues_u, cellvalues_ϕ, cellvalues_θ, qe, material, state, state_old)
    nbase_u = getnbasefunctions(cellvalues_u)
    nbase_ϕ = getnbasefunctions(cellvalues_ϕ)
    nbase_θ = getnbasefunctions(cellvalues_θ)
    u▄, ϕ▄, θ▄ = 1, 2, 3
    reinit!(cellvalues_u, cell)
    reinit!(cellvalues_ϕ, cell)
    reinit!(cellvalues_θ, cell)
    ϕe = qe[nbase_u + 1 : nbase_u + nbase_ϕ]
    θe = qe[nbase_u + nbase_ϕ + 1 : end]
    ue = qe[1 : nbase_u] 

    Gc = material.Gc
    ℓ = material.ℓ

    for q_point in 1:getnquadpoints(cellvalues_u)
        
        ϕ = state[q_point].ϕ;
        θ = state[q_point].θ;
        # ϕ = function_value(cellvalues_ϕ, q_point, ϕe);
        # θ = function_value(cellvalues_θ, q_point, θe);
        ∇θ = function_gradient(cellvalues_θ, q_point, θe);

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
            for j in 1:nbase_θ
                θ′ = shape_value(cellvalues_θ, q_point, j);
                Ke[BlockIndex((u▄, θ▄), (i, j))] -= (δεu ⊡ Dᵉ ⊡ one(SymmetricTensor{2, 3})* material.α * θ′) * dΩᵤ
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
                if state[q_point].H >= state_old[q_point].H
                    Ke[BlockIndex((ϕ▄, u▄), (i, j))] += (gd′ *δϕ * σ⁺ ⊡ ε̄u) * dΩᵤ
                end   
            end
            for j in 1:nbase_θ
                θ′ = shape_value(cellvalues_θ, q_point, j);
                if state[q_point].H > state_old[q_point].H
                    Ke[BlockIndex((ϕ▄, θ▄), (i, j))] -= (gd′ *δϕ * σ⁺ ⊡ one(SymmetricTensor{2, 3})* material.α * θ′) * dΩᵤ
                end   
            end
        end
        for i in 1:nbase_θ
            δθ = shape_value(cellvalues_θ, q_point, i);
            ∇δθ = shape_gradient(cellvalues_θ, q_point, i);
            for j in 1:nbase_θ
                θ′ = shape_value(cellvalues_θ, q_point, j);
                ∇θ′ = shape_gradient(cellvalues_θ, q_point, j);
                Ke[BlockIndex((θ▄, θ▄), (i, j))] += (ρ * Cp/Δt * δθ * θ′  + gd * k₀ * ∇δθ ⋅ ∇θ′) * dΩᵤ
            end
            # for j in 1:nbase_u
            #     Ke[BlockIndex((θ▄, u▄), (i, j))] += 0.0
            # end
            for j in 1:nbase_ϕ
                ϕ′ = shape_value(cellvalues_ϕ, q_point, j)
                ∇ϕ′ = shape_gradient(cellvalues_ϕ, q_point, j);
                Ke[BlockIndex((θ▄, ϕ▄), (i, j))] += (gd′ * k₀* ∇δθ ⋅ ∇θ *ϕ′) * dΩᵤ
            end
        end
    end    
end

function assemble_global_K(q, K, cellvalues_u,cellvalues_ϕ, cellvalues_θ, dh, material, states, states_old)

    assembler = start_assemble(K)
    nu = getnbasefunctions(cellvalues_u)
    nϕ = getnbasefunctions(cellvalues_ϕ)
    nθ = getnbasefunctions(cellvalues_θ)

    Ke = BlockArray(zeros(nu + nϕ + nθ, nu + nϕ + nθ), [nu, nϕ, nθ], [nu, nϕ, nθ]) # local stiffness matrix
    for (i, cell)  in enumerate(CellIterator(dh))
        fill!(Ke, 0)
        eldofs = celldofs(cell)
        qe = q[eldofs]
        state = @view states[:, i]
        state_old = @view states_old[:, i]
        assemble_element_ke!(Ke, cell, cellvalues_u, cellvalues_ϕ, cellvalues_θ, qe, material, state, state_old)

        assemble!(assembler, eldofs, Ke)
    end
    return K
end

# function CalculateRho(q, K, cellvalues_u,cellvalues_ϕ, cellvalues_θ, dh, material, states, states_old)
#     undofs = dh_u.ndofs;
#     ϕndofs = dh_ϕ.ndofs;
#     K = assemble_global_K(q, K, cellvalues_u,cellvalues_ϕ, cellvalues_θ, dh, material, states, states_old)
#     A = copy(K)
#     A[undofs+1 : ϕndofs+undofs,1:undofs] .=0.0;
#     B = allocate_matrix(dh)
#     B[1:undofs, undofs+1 : ϕndofs+undofs] .= @view K[1:undofs,1:ϕndofs]
    
#     return A, B
# end