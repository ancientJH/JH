function symmetrize_lower!(K)
    for i in 1:size(K,1)
        for j in i+1:size(K,1)
            K[i,j] = K[j,i]
        end
    end
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
        else
            assemble_element_θ!(Ke, Re, cellvalues, qe, material, state, state_old);
        end
        assemble!(assembler, eldofs, Ke, Re);
        
    end
    return K, R
end;

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
        else
            assemble_residual_θ!(Re, cellvalues, qe, material, state, state_old, store);
        end
        R[eldofs] += Re;
    end
    return R
end;