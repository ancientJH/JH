
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

function assemble_global_ru_withϕ(u::Vector, ϕ::Vector, cellvalues_u, cellvalues_ϕ, dh, dh_u::DofHandler, dh_ϕ::DofHandler, material::Material, states, states_old, tag::String, store::Bool)
    nbase = getnbasefunctions(cellvalues_u);
    Re = zeros(nbase);
    R = zeros(ndofs(dh_u));
    @inbounds for (i, cell) in enumerate(CellIterator(dh))
        reinit!(cellvalues_u, cell);
        reinit!(cellvalues_ϕ, cell);
        fill!(Re, 0);
        eldofs_u_local = celldofs(dh_u, i)
        eldofs_ϕ_local = celldofs(dh_ϕ, i)
        ue = u[eldofs_u_local];
        ϕe = ϕ[eldofs_ϕ_local];
        state = @view states[:, i];
        state_old = @view states_old[:, i];
        assemble_residual_u_withϕ!(Re, cellvalues_u, cellvalues_ϕ, ue, ϕe, material, state, state_old, store);
        # R[eldofs_u_local] += Re;
        # 使用 addindex! 将局部残差 Re 添加到全局残差 R
        @inbounds for (local_idx, global_idx) in pairs(eldofs_u_local)
            Ferrite.addindex!(R, Re[local_idx], global_idx)
        end
    end
    return R
end;

function assemble_global_rϕ_withu(u::Vector, ϕ::Vector, cellvalues_u, cellvalues_ϕ, dh, dh_u::DofHandler, dh_ϕ::DofHandler, material::Material, states, states_old, tag::String, store::Bool)
    nbase = getnbasefunctions(cellvalues_ϕ);
    Re = zeros(nbase);
    R = zeros(ndofs(dh_ϕ));

    @inbounds for (i, cell) in enumerate(CellIterator(dh))
        reinit!(cellvalues_u, cell);
        reinit!(cellvalues_ϕ, cell);
        fill!(Re, 0);
        eldofs_u_local = celldofs(dh_u, i)
        eldofs_ϕ_local = celldofs(dh_ϕ, i)
        ue = u[eldofs_u_local];
        ϕe = ϕ[eldofs_ϕ_local];
        state = @view states[:, i];
        state_old = @view states_old[:, i];
        assemble_residual_ϕ_withu!(Re, cellvalues_u, cellvalues_ϕ, ue, ϕe, material, state, state_old, store);
        # R[eldofs_ϕ_local] += Re;
        # 使用 addindex! 将局部残差 Re 添加到全局残差 R
        @inbounds for (local_idx, global_idx) in pairs(eldofs_ϕ_local)
            Ferrite.addindex!(R, Re[local_idx], global_idx)
        end
    end
    return R
end;
