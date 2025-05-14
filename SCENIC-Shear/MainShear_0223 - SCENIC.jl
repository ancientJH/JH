
##########################################################
using BlockArrays
using Ferrite
using FerriteGmsh
using HDF5
using IterativeSolvers
using JLD2
using LinearAlgebra
using MAT
using SparseArrays
using Tensors
using WriteVTK
using Printf
using Arpack
using Profile
using ProfileView
# using SuiteSparse
# using Metis
#using infiltrator
#Infiltrator.toggle_async_check(false)  # 禁用异步检查
#Infiltrator.clear_disabled!() 启用异步检查
@enum StrainDecomp Isotropic VolDev Spectral
@enum DegradType QuadraticDegradation WuDegradation HughesDegradation
include("Subfunction - SCENIC.jl")
##########################################################
function Problem(solver, output, cellvalues_ϕ, cellvalues_u,
	dh, dh_u, dh_ϕ, bch, ch_u, ch_ϕ,grid, Mat)
	rho_ref = 0.25
	
	u = zeros(ndofs(dh_u))
	ϕ = zeros(ndofs(dh_ϕ))

	K_u = create_sparsity_pattern(dh_u)
	K_ϕ = create_sparsity_pattern(dh_ϕ)
	K = create_sparsity_pattern(dh)

	nqp = getnquadpoints(cellvalues_u)
	states = [HistoryVariable() for _ in 1:nqp, _ in 1:getncells(grid)] #初始化
	states_old = [HistoryVariable() for _ in 1:nqp, _ in 1:getncells(grid)] # 初始化
	Disp = 0

	output_folder = "Shear20390ele"
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
	
	for timestep in start_timestep:length(loadsteps)  #timestep in 0:n_timesteps
		solve_time = @elapsed begin
		flag = "Staggered"
		flag_eta = false  # 初始化 flag_eta
		Disp = loadsteps[timestep]       #Disp = timestep * solver.u_max / n_timesteps
		update!(ch_u, Disp)
		update!(ch_ϕ)
		#update!(bch, Disp)
		apply!(u, ch_u)
		apply!(ϕ, ch_ϕ)
		#apply!(q, bch)
        iterations_out = 0
		Ru_first = 1.0
		Rϕ_first = 1.0
		recent_norms = Float64[]
		output.totalIterations_u = 0; output.totalIterations_ϕ = 0; output.totalIterations_outer = 0;
		for newton_itr ∈ 1:solver.nitr_outer
            iterations_out += 1
			if newton_itr > solver.nitr_outer
				error("Reached maximum Newton iterations, aborting")
				break
			end

			#_, r_ϕ = assemble_global!(u, ϕ, K, cellvalues_u, cellvalues_ϕ, dh, dh_u, dh_ϕ, Mat, states, states_old, "ϕ");
			r_u = assemble_global_r(u, cellvalues_u, dh_u, Mat, states, states_old, "u");
            r_ϕ = assemble_global_r(ϕ, cellvalues_ϕ, dh_ϕ, Mat, states, states_old, "ϕ");
			R_u = r_u[Ferrite.free_dofs(ch_u)]
			R_ϕ = r_ϕ[Ferrite.free_dofs(ch_ϕ)]
			norm_rϕ = norm(R_ϕ) #norm_r = maximum(abs.(r_ϕ[Ferrite.free_dofs(ch_ϕ)]))
			norm_ru = norm(R_u)

			push!(recent_norms, norm_ru)
            if length(recent_norms) > 3
                popfirst!(recent_norms)  # 保持数组长度为 3
            end
			if length(recent_norms) > 1
				if flag == "Staggered" && (recent_norms[end]/recent_norms[end - 1] > rho_ref) 
					flag = "SCENIC"
					recent_norms = Float64[]
					flag_eta = false;
				end
			end

			if length(recent_norms) > 1
				if flag == "SCENIC" && (recent_norms[end] > recent_norms[end - 1]) 
					flag_eta = true;
				end
				
			end
			
            if newton_itr == 1
				Ru_first = norm_ru
				Rϕ_first = norm_rϕ
			end
			println("norm_ru=$norm_ru,norm_rϕ=$norm_rϕ,Ru_first=$Ru_first,Rϕ_first=$Rϕ_first ");
			if norm_rϕ <= solver.TOL_ϕ*Rϕ_first && norm_ru <= solver.TOL_u*Ru_first
				print("\n Disp = $Disp, converged in $(newton_itr-1) outer iterations\n")
				save_results(output_folder,loadsteps[timestep], u, ϕ, states, states_old)	
				output.totalIterations_outer = newton_itr - 1	
				states_old .= states
				break
			end
			
			print("\n REBUILDING STIFFNESS at Disp = $Disp, current iteration count is $(newton_itr-1)\n")

			if flag == "SCENIC"
				# Profile.clear()
				Δu =  Masterstep_u!(u, ϕ, K, R_u, R_ϕ, cellvalues_u, cellvalues_ϕ, dh, dh_u, dh_ϕ, bch, ch_u, ch_ϕ, Mat, states, states_old, flag_eta)
			    output.totalIterations_u += 1 
				# ProfileView.view()
				# Profile.print() 
				# error("stop")
                u += Δu
				states = updatePsi(cellvalues_u, u, dh_u, Mat, states, states_old); # 更新一下 H
			elseif flag == "Staggered"
				u, nitr_u = Newton_raphson!(u, K_u, cellvalues_u, dh_u, ch_u, grid, Mat, states, states_old, Ru_first, Rϕ_first, "u");
			    output.totalIterations_u += nitr_u
			end
            

			ϕ, nitr_ϕ = Newton_raphson!(ϕ, K_ϕ, cellvalues_ϕ, dh_ϕ, ch_ϕ, grid, Mat, states, states_old, Ru_first, Rϕ_first, "ϕ");
			output.totalIterations_ϕ += nitr_ϕ
			
		end
		end #时间的end
		if timestep in (0:output.plotFrequency:length(loadsteps)) || timestep == 1 || timestep == length(loadsteps)
			timestep_str = replace(@sprintf("%.5f", loadsteps[timestep]), "." => "_")
			simtime = timestep
			saved_file = paraview_collection(joinpath(output_folder, "TimeSeries"); append = true) do pvd
				vtk_grid(joinpath(output_folder, "time_$timestep_str"), dh_ϕ) do vtk
					vtk_point_data(vtk, dh_ϕ, ϕ)
					vtk_point_data(vtk, dh_u, u)
					# vtk_point_data(vtk, σ_projected, "σ(MPa)")
					pvd[simtime] = vtk
				end
			end
			output.plotframe += 1
		end
		#Write history output
		if timestep in (0:output.historyFrequency:length(loadsteps)) || timestep == 1 || timestep == length(loadsteps)
			# 计算输出力
			F_x, F_y = OutputForce(u, cellvalues_u, dh_u, grid, Mat, states, states_old, output.OutputSet, "u")
			# 确保输出文件夹存在
			history_folder = output_folder  # 指定输出文件夹
			if !isdir(history_folder)
				mkdir(history_folder)
			end
			# 构造历史输出文件路径
			history_file_path = joinpath(history_folder, "HistoryOutput.txt")
			# 打开文件并写入历史数据
			historyfile = open(history_file_path, "a")
			write(historyfile, "$timestep, $(output.totalIterations_u), $(output.totalIterations_ϕ),  " *
							   "$(output.totalIterations_outer),$solve_time, $Disp, $F_x, $F_y\n")
			close(historyfile)
		end
		# states_old .= states
	end
	print("Computation complete! Runtime information:\n")
end;
##########################################################
dim = 2;
ElementOrder = 1;
QuadratureOrder = 2;
#定义插值形函数
# ElementShape = RefTriangle;
# ip = Lagrange{ElementShape, ElementOrder}()
# qr = QuadratureRule{ElementShape}(QuadratureOrder)
# #Add machinery for cell (element) values
# cellvalues_u = CellValues(qr, ip^2)
# cellvalues_ϕ = CellValues(qr, ip)
ElementShape = RefTetrahedron;
ip = Lagrange{dim, ElementShape, ElementOrder}()
qr = QuadratureRule{dim, ElementShape}(QuadratureOrder)

cellvalues_u = CellVectorValues(qr, ip)
cellvalues_ϕ = CellScalarValues(qr, ip)
##########################################################
grid = togrid("Shear20390ele.geo");
addnodeset!(grid, "top", x -> x[2] ≈ 1.0);
# dbc₁ = Dirichlet(:u, getfacetset(grid, "top"), (x, t) -> t, 2);
# dbc₂ = Dirichlet(:u, getfacetset(grid, "bottom"), (x, t) -> [0, 0], [1, 2]);
# dbc₁₂ = Dirichlet(:u, getnodeset(grid, "CornerPoint"), (x, t) -> [0, 0], [1, 2]);
# dbc₃ = Dirichlet(:ϕ, getnodeset(grid, "Crack"), (x, t) -> 1); 
dbc₁ = Dirichlet(:u, getfaceset(grid, "top"), (x, t) -> [t,0], [1, 2]);
dbc₂ = Dirichlet(:u, getfaceset(grid, "bottom"), (x, t) -> [0, 0], [1, 2]);

dh = MixedDofHandler(grid);
#dh = DofHandler(grid);
add!(dh, :u, dim);
add!(dh, :ϕ, 1);
close!(dh);
renumber!(dh, DofOrder.FieldWise()); #这句会让自由度的编号为：u所有的自由度，d所有的自由度（针对单元内部）

bch = ConstraintHandler(dh);
add!(bch, dbc₁);
add!(bch, dbc₂); #add!(ch_u, dbc₁₂); 
close!(bch);
update!(bch,0.0)

dh_u = DofHandler(grid);
add!(dh_u, :u, dim);
close!(dh_u);

ch_u = ConstraintHandler(dh_u);
add!(ch_u, dbc₁);
add!(ch_u, dbc₂); #add!(ch_u, dbc₁₂); 
close!(ch_u);
update!(ch_u, 0.0);

dh_ϕ = DofHandler(grid);
add!(dh_ϕ, :ϕ, 1);
close!(dh_ϕ);

ch_ϕ = ConstraintHandler(dh_ϕ); #add!(ch_ϕ, dbc₃);
close!(ch_ϕ);
update!(ch_ϕ, 0.0);

##########################################################
loadsteps = vcat(collect(0.0005:0.0005:0.009), collect(0.0092:0.0002:0.03))
nitr_inner = 2000;  #牛顿迭代最大次数
nitr_outer = 2000; #外部迭代次数
LRuactive = length(free_dofs(ch_u));
LRϕactive = length(free_dofs(ch_ϕ));
TOL_u = 1e-5;
TOL_ϕ = 1e-5;
# TOL_Ru = sqrt(LRuactive)*1e-6;
# TOL_Rϕ = sqrt(LRϕactive)*1e-6;
##########################################################
plot_frequency = 1;
history_output_frequency = 1;
a₀ = 0.0;
CrackDir = 1;
outputset = "top";
##########################################################
E = 210e3;#MPa
ν = 0.3;
Gc = 2.7; #N/mm=kJ/m^2
σc = 2445.42; #MPa
ℓ₀ = 0.01;
s = 200.0;
##########################################################
mat = Material(E, ν, Gc, σc, ℓ₀, s, Spectral, QuadraticDegradation, dim);
solver = create_solver_state(loadsteps, nitr_inner, nitr_outer, TOL_u, TOL_ϕ);
output = OutputVariables(plot_frequency, history_output_frequency, a₀, CrackDir, outputset);
@time Problem(solver, output, cellvalues_ϕ, cellvalues_u, dh, dh_u, dh_ϕ, bch, ch_u, ch_ϕ, grid, mat);