function Newton_raphson!(q::Vector, K::SparseMatrixCSC, cellvalues, dh::DofHandler, ch::ConstraintHandler, grid::Grid, Mat::Material, states, states_old, Ru_first, Rϕ_first, Rθ_first, tag::String)
    iterations = 0;
    for nitr = 1:(solver.nitr_inner+1);
        if nitr > solver.nitr_inner;
            error("Reached maximum Newton iterations, aborting");
            break;
        end;
        K, r = assemble_global(q, cellvalues, K, dh, Mat, states, states_old, tag::String);
        apply_zero!(K, r, ch);
        norm_r = norm(r[Ferrite.free_dofs(ch)]);  # norm_r = maximum(abs.(r[Ferrite.free_dofs(ch)]))
        if tag =="u"
            TOL = solver.TOL_u * Ru_first;
        elseif tag =="ϕ"
            TOL = solver.TOL_ϕ * Rϕ_first;
        else
            TOL = solver.TOL_θ * Rθ_first;
        end
        if (norm_r < TOL) && (nitr > 1);
            break;
        end
        iterations += 1;
        K_active = K[Ferrite.free_dofs(ch),Ferrite.free_dofs(ch)];
		r_active =  r[Ferrite.free_dofs(ch)];
        # Δq = K_active \ r_active;
        q[Ferrite.free_dofs(ch)] -= K_active \ r_active;
    end
        print(tag*" converged in $iterations iterations \n");
    return q, iterations, K
end
