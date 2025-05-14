using Ferrite, Tensors, LinearAlgebra, BlockArrays, SparseArrays, FerriteGmsh, WriteVTK, Printf, Arpack
using MAT, JLD, HDF5
@enum StrainDecomp Isotropic VolDev Spectral
@enum DegradType QuadraticDegradation WuDegradation HughesDegradation

include("../Src/Material_Hisvar_Solver_Output.jl")
include("../Src/Assemble/Assemble_u.jl")
include("../Src/Assemble/Assemble_d.jl")
include("../Src/Assemble/Assemble_T.jl")
include("../Src/Assemble/Assemble_Global.jl")
include("../Src/Assemble/Assemble_DOC.jl")
include("../Src/Assemble/Assemble_ud.jl")
include("../Src/Newton_Force_Crack.jl")
include("../Src/Constitutive_Degradation.jl")
include("../Src/SCENIC.jl")
function Problem(solver, output, cellvalues_ϕ, cellvalues_θ, cellvalues_u, dh_ϕ, dh_θ, dh_u, dh, ch_ϕ, ch_θ, ch_u, bch, grid, Mat)
    
    u = zeros(ndofs(dh_u));
    ϕ = zeros(ndofs(dh_ϕ));
    θ = Mat.θ₀ * ones(ndofs(dh_θ));
    K_u = create_sparsity_pattern(dh_u);
    K_ϕ = create_sparsity_pattern(dh_ϕ);
    K_θ = create_sparsity_pattern(dh_θ);
    K = create_sparsity_pattern(dh);

    nqp = getnquadpoints(cellvalues_u);
    states = [HistoryVariable() for _ in 1:nqp, _ in 1:getncells(grid)];
    states_old = [HistoryVariable() for _ in 1:nqp, _ in 1:getncells(grid)];

    output_folder = "Thermal2_3wele"
	if !isdir(output_folder)
		mkdir(output_folder)
	end
	start_timestep = 1
    if isfile(joinpath(output_folder, "time_0_01060.jld"))  # 检查文件是否存在
		timestep, u, ϕ, states, states_old = load_results(output_folder, 0.01060)  # 从文件夹 "Shear" 加载时间步 0.001
		start_timestep = findfirst(x -> x == timestep, loadsteps)  # 找到对应的时间步索引
		println("Resuming from timestep $timestep (index $start_timestep)")
		start_timestep = start_timestep + 1  # 从下一个时间步开始
	else
		println("Starting from the beginning")
	end

    for timestep in 1: length(loadsteps)
        solve_time = @elapsed begin
        rho_ref1 = 0.1;
        rho_ref2 = 0.1;
        flag1 = "Staggered";
        flag2 = "Staggered";
        flag_eta1 = false;
        flag_eta2 = false;
        recent_norms_θ = Float64[]
		recent_norms_u = Float64[]
        apply!(u, ch_u);
        apply!(ϕ, ch_ϕ);
        apply!(θ, ch_θ);
        update!(bch, θ₀);
        apply!(vcat(u, ϕ, θ), bch);
        if timestep==1
            for (i, cell) in enumerate(CellIterator(dh_θ))
                reinit!(cellvalues_θ, cell);
                eldofs = celldofs(cell);
                θe = θ[eldofs];
                for q_point in 1:nqp
                    θᵢ = function_value(cellvalues_θ, q_point, θe);
                    states[q_point, i] = HistoryVariable(0.0, 0.0, θᵢ);
                    states_old[q_point, i] = HistoryVariable(0.0, 0.0, θᵢ);
                end
            end
        end
        iterations_out = 0
        iterations_inter = 0
	    Ru_first = 1.0
	    Rϕ_first = 1.0
        Rθ_first = 1.0
        for master_iter1 in 1:solver.nitr_outer
            iterations_out +=1
            if master_iter1 > solver.nitr_outer
                error("Reached maximum Newton iterations, aborting")
                break
            end
            # states = updatePsi(cellvalues_u, u, dh_u, Mat, states, states_old);
            r_θ = assemble_global_r(θ, cellvalues_θ, dh_θ, Mat, states, states_old, "θ");
            r_u = assemble_global_r(u, cellvalues_u, dh_u, Mat, states, states_old, "u");
            r_ϕ = assemble_global_r(ϕ, cellvalues_ϕ, dh_ϕ, Mat, states, states_old, "ϕ");
            

            norm_rϕ = norm(r_ϕ[free_dofs(ch_ϕ)]);
            norm_ru = norm(r_u[free_dofs(ch_u)]);
            norm_rθ = norm(r_θ[free_dofs(ch_θ)]);

            if master_iter1==1
                Ru_first = max(norm_ru,1.0)
	            Rϕ_first = max(norm_rϕ,1e-4)
                Rθ_first = max(norm_rθ,1e-4)
            end

            push!(recent_norms_θ, norm_rθ)
            if length(recent_norms_θ) > 3
                popfirst!(recent_norms_θ)  # 保持数组长度为 3
            end
			if length(recent_norms_θ) > 1
				if flag1 == "Staggered" && (recent_norms_θ[end]/recent_norms_θ[end - 1] > rho_ref1) 
					flag1 = "SCENIC"
					recent_norms_θ = Float64[]
					flag_eta1 = false;
				end
			end

            if length(recent_norms_θ) > 1
				if flag1 == "SCENIC" && (recent_norms_θ[end] > recent_norms_θ[end - 1]) 
					flag_eta1 = true;
				end
			end

            print("Time step @time = $timestep, master_iter1 = $master_iter1, Ru = $norm_ru, Rϕ = $norm_rϕ, Rθ = $norm_rθ\n");
            if norm_rϕ <= solver.TOL_ϕ * Rϕ_first &&norm_ru <= solver.TOL_u * Ru_first &&norm_rθ <= solver.TOL_θ * Rθ_first
                print("\n Time step @time = $timestep, Ru, Rϕ and RT converged 1st in $(master_iter1) master iterations\n");
                break;
            end
            
            
            ######################  update θ fields  #############################
            if flag1 == "Staggered"
                θ, nitr_θ, K_θ = Newton_raphson!(θ, K_θ, cellvalues_θ, dh_θ, ch_θ, grid, Mat, states, states_old, Ru_first, Rϕ_first, Rθ_first, "θ");
                output.totalIterations_θ += nitr_θ; # 不用更新 H 了，因为后面计算 r_u 的时候会更新 H
            elseif flag1 == "SCENIC"
                Δθ = Masterstep_θ!(u, ϕ, θ, K, r_u, r_ϕ, r_θ, cellvalues_u, cellvalues_ϕ, cellvalues_θ, dh, dh_u, dh_ϕ, dh_θ, bch, ch_u, ch_ϕ, ch_θ, Mat, states, states_old, flag_eta1)
			    output.totalIterations_θ += 1 # 无法避免计算本构
                θ += Δθ # 这里要将 Δθ 加到 states 里，以及更新 H
                println("θ have been updated by SCENIC")
                states = updateT(cellvalues_θ, θ, dh_θ, Mat, states, states_old); 
                states = updatePsi(cellvalues_u, u, dh_u, Mat, states, states_old); # 更新一下 H
            end

            for master_iter2 in 1:solver.nitr_outer
                iterations_inter +=1
                if master_iter2 > solver.nitr_outer
                    error("Reached maximum Newton iterations, aborting")
                    break
                end

                r_u = assemble_global_r(u, cellvalues_u, dh_u, Mat, states, states_old, "u"); # 要先算 r_u 再算 r_θ 因为 r_θ 需要用到 H
                r_ϕ = assemble_global_r(ϕ, cellvalues_ϕ, dh_ϕ, Mat, states, states_old, "ϕ");
                

                norm_ru = norm(r_u[free_dofs(ch_u)]);
                norm_rϕ = norm(r_ϕ[free_dofs(ch_ϕ)]);

                push!(recent_norms_u, norm_ru)
				if length(recent_norms_u) > 3
					popfirst!(recent_norms_u)  # 保持数组长度为 3
				end
				if length(recent_norms_u) > 1
					if flag2 == "Staggered" && (recent_norms_u[end]/recent_norms_u[end - 1] > rho_ref2) 
						flag2 = "SCENIC"
						recent_norms_u = Float64[]
						flag_eta2 = false;
					end
				end

				if length(recent_norms_u) > 1
					if flag2 == "SCENIC" && (recent_norms_u[end] > recent_norms_u[end - 1]) 
						flag_eta2 = true;
					end
				end

                print("Time step @time = $timestep, master_iter2 = $master_iter2, Ru = $norm_ru, Rϕ = $norm_rϕ\n");
                if norm_ru <= solver.TOL_u * Ru_first &&norm_rϕ <= solver.TOL_ϕ * Rϕ_first
                    output.totalIterations_inter += iterations_inter;
                    print("\n Time step at master_iter1 = $master_iter1, Ru and Rϕ converged in $(master_iter2) inner iterations\n");
                    break;
                end
                ######################  update u fields  #############################
                if flag2 == "Staggered"
                    u, nitr_u, K_u = Newton_raphson!(u, K_u, cellvalues_u, dh_u, ch_u, grid, Mat, states, states_old, Ru_first, Rϕ_first, Rθ_first, "u");
                    output.totalIterations_u += nitr_u; 
                elseif flag2 == "SCENIC"
                    Δu = Masterstep_u!(u, ϕ, θ, K, r_u, r_ϕ, r_θ, cellvalues_u, cellvalues_ϕ, cellvalues_θ, dhud, dh_u, dh_ϕ, dh_θ, bchud, ch_u, ch_ϕ, ch_θ, Mat, states, states_old, flag_eta1)
                    output.totalIterations_u += 1 # 无法避免计算本构
                    u += Δu
                    println("u have been updated by SCENIC")
                    states = updatePsi(cellvalues_u, u, dh_u, Mat, states, states_old); # 更新一下 H
                end
                
                ######################  update ϕ fields  #############################
                ϕ, nitr_ϕ, K_ϕ = Newton_raphson!(ϕ, K_ϕ, cellvalues_ϕ, dh_ϕ, ch_ϕ, grid, Mat, states, states_old, Ru_first, Rϕ_first, Rθ_first, "ϕ");
                output.totalIterations_ϕ += nitr_ϕ;
                
            end
            output.totalIterations_outer = iterations_out
        end
        end #时间的end
        if timestep in (0:output.plotFrequency:length(loadsteps)) || timestep == 1 || timestep == length(loadsteps)
            
            simtime = timestep;
            timestep_str = replace(@sprintf("%.5f", timestep), "." => "_")
            saved_file = paraview_collection(joinpath(output_folder, "TimeSeries"); append = true) do pvd
				vtk_grid(joinpath(output_folder, "time_$timestep_str"), dh_ϕ) do vtk
					vtk_point_data(vtk, dh_ϕ, ϕ)
					vtk_point_data(vtk, dh_u, u)
                    vtk_point_data(vtk, dh_θ, (θ .-273))
					pvd[simtime] = vtk
				end
			end
            K = assemble_global_K(vcat(u, ϕ, θ), K, cellvalues_u,cellvalues_ϕ, cellvalues_θ, dh, Mat, states, states_old)
            save("Examples2/T_ud_0224_Testing/Data_$(timestep).jld", "u", u, "ϕ", ϕ, "θ", θ, "states", states, "states_old", states_old,
            "K_θ",K_θ,"K_ϕ",K_ϕ,"K_u",K_u, "K",K)
            output.plotframe += 1
        end
        if timestep in (0:output.historyFrequency:length(loadsteps)) || timestep == 1 || timestep == length(loadsteps);
            history_folder = output_folder  # 指定输出文件夹
			if !isdir(history_folder)
				mkdir(history_folder)
			end
			# 构造历史输出文件路径
			history_file_path = joinpath(history_folder, "HistoryOutput.txt")
            historyfile = open(history_file_path, "a")
            write(historyfile, "$timestep, $(output.totalIterations_u), $(output.totalIterations_ϕ), $(output.totalIterations_θ),  " *
                               "$(output.totalIterations_outer), $(output.totalIterations_inter), $solve_time\n")
            close(historyfile)
        end
        states_old .= states;
    end
    print("Computation complete! Runtime information:\n")
end

##########################################################
if isdir("Examples2/T_ud_0224_Testing")
    rm("Examples2/T_ud_0224_Testing",recursive=true)     
else
    mkpath("Examples2/T_ud_0224_Testing")
end
##########################################################
E = 370e9;#Pa
ν = 0.3;
Gc= 42.47;  #N/m=J/m^2
σc= 180e6;  #Pa
ℓ₀ = 7.5e-5;  #m
s = 200.; 
ρ = 3980.;  # kg/m^3
Cp = 880.;   # J/(kg⋅ K) specific heat capacity   
k₀ = 31.;     # W/(m⋅ K) thermal conductivity  W=J/s
α = 7.5e-6; # 1/K thermal expansion coefficient
θ₀ = 573.;
Δt = 1e-4;
##########################################################
dim = 2;
ElementOrder = 1;
QuadratureOrder = 2;
ElementShape = RefTetrahedron;
ip = Lagrange{dim, ElementShape, ElementOrder}()
qr = QuadratureRule{dim, ElementShape}(QuadratureOrder)
cellvalues_u = CellVectorValues(qr, ip)
cellvalues_ϕ = CellScalarValues(qr, ip)
cellvalues_θ = CellScalarValues(qr, ip)
##########################################################
grid = togrid("../MeshFiles/Thermal.geo");
addnodeset!(grid, "top", x -> x[2] ≈ 5e-3);
dbc_u1 = Dirichlet(:u, getfaceset(grid, "top"),   (x, t) -> 0, 2);  
dbc_u2 = Dirichlet(:u, getfaceset(grid, "right"), (x, t) -> 0, 1);
dbc_θ1 = Dirichlet(:θ, getfaceset(grid, "left"),  (x, t) ->293); 
dbc_θ2 = Dirichlet(:θ, getfaceset(grid, "bottom"),(x, t) ->293); 
##########################################################
dh_u = DofHandler(grid); add!(dh_u, :u, dim);  close!(dh_u);
ch_u = ConstraintHandler(dh_u); 
add!(ch_u, dbc_u1); add!(ch_u, dbc_u2); 
close!(ch_u); update!(ch_u, 0.0);

dh_ϕ = DofHandler(grid); add!(dh_ϕ, :ϕ, 1); close!(dh_ϕ);
ch_ϕ = ConstraintHandler(dh_ϕ);
close!(ch_ϕ); update!(ch_ϕ, 0.0);

dh_θ = DofHandler(grid); add!(dh_θ, :θ, 1); close!(dh_θ);
ch_θ = ConstraintHandler(dh_θ);
add!(ch_θ, dbc_θ1); add!(ch_θ, dbc_θ2); 
close!(ch_θ); update!(ch_θ, θ₀);
##########################################################
dh = MixedDofHandler(grid);add!(dh,:u, dim); add!(dh, :ϕ, 1); add!(dh, :θ, 1); close!(dh);
renumber!(dh, DofOrder.FieldWise());
dhud = MixedDofHandler(grid);add!(dhud,:u, dim); add!(dhud, :ϕ, 1);  close!(dhud);
renumber!(dhud, DofOrder.FieldWise());
bch = ConstraintHandler(dh); add!(bch, dbc_u1); add!(bch, dbc_u2); add!(bch, dbc_θ1); add!(bch, dbc_θ2); close!(bch); 
bchud = ConstraintHandler(dhud); add!(bchud, dbc_u1); add!(bchud, dbc_u2); close!(bchud); 
##########################################################
loadsteps = collect(LinRange(1, 100, 100));   # collect(1:1:50);
# loadsteps = collect(1:1:50);
nitr_inner = 2000;  #牛顿迭代最大次数
nitr_outer = 2000; #外部迭代次数
TOL_Ru = 1e-6;
TOL_Rϕ = 1e-6;
TOL_Rθ = 1e-6;
##########################################################
plot_frequency = 1;
history_output_frequency = 1;
a₀ = 0.0;
CrackDir = 1;
outputset = "top";
##########################################################
mat = Material(E, ν, Gc, σc, ℓ₀, ρ, Cp, k₀, α, θ₀, Δt, s, Spectral , QuadraticDegradation, dim);
solver = SolverState(loadsteps, nitr_inner, nitr_outer, TOL_Ru, TOL_Rϕ, TOL_Rθ);
output = OutputVariables(plot_frequency, history_output_frequency, a₀, CrackDir, outputset);
@time Problem(solver, output, cellvalues_ϕ, cellvalues_θ, cellvalues_u, dh_ϕ, dh_θ, dh_u, dh, ch_ϕ, ch_θ, ch_u, bch, grid, mat);
