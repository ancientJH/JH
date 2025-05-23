
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
        end
        assemble!(assembler, eldofs, Ke, Re);
        
    end
    return K, R
end;