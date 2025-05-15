%% 准备工作
close all;
clear all;
clc;
outputfile = 'Shear25422ele'; %需要修改这里，网格文件需要也是这个名字
if exist(outputfile, 'dir') ~= 0
    error('Directory already exists: %s', fullfile(outputfile, '*.*')); 
else
    mkdir(outputfile); 
end
%% 问题定义
load_steps = [0.5:0.5:9,9.2:0.2:30]; % 定义加载量
elementType = 'P12D';problemType = 'PhaseField'; body_force = @user_body_force; traction = @user_traction;
domain = [0,1; 0,1]; nDim = 2; % 三角形网格，二维
nDoF = nDim+1; DirichletBCs = @PhaseFieldDirichletBCs; constitutive = @ModelC; % 模型和边界条件
load([outputfile '.mat'], 's', 't', 'p'); 
nElements = size(t,2);nNodes = size(p,2);nNodesElement = 3; Coord = p;
IEN = t(1:nNodesElement, 1:nElements); nEquations = nNodes * nDoF;
base = (IEN - 1) * nDoF; LM = reshape(bsxfun(@plus, reshape(base, 1, []), (1:nDoF)'), [], nElements);
%% 可视化网格
% vertices = Coord'; faces = IEN'; fig = figure(1); fig.Renderer = 'OpenGL'; % 使用硬件加速渲染
% p = patch('Faces', faces,'Vertices', vertices,'FaceVertexCData', ones(size(faces,1),1),... % 每个单元一个颜色值
%     'FaceColor', 'flat','EdgeColor', 'k','LineWidth', 0.5);     
% axis equal tight; colormap(parula); c = colorbar('southoutside'); c.Label.String = 'Element Data';
% set(gca, 'Visible', 'off');set(gca, 'SortMethod', 'depth'); drawnow limitrate; 
%% 边界条件
delta_u = 1e-3; tolerance = 1e-6; nDofPerNode = size(Coord, 1) + 1; nNodes = size(Coord, 2);
BCIndices = false(nDofPerNode, nNodes);BCVal = zeros(nDofPerNode, nNodes, 'single');  
FIndices = zeros(nDofPerNode-1, nNodes, 'logical'); y_coord = Coord(2,:); 
bottom_mask = abs(y_coord) <= tolerance;BCIndices(1:2, bottom_mask) = true;BCVal(1:2, bottom_mask) = 0;
top_mask = abs(y_coord - 1) <= tolerance; BCIndices(2, top_mask) = true;BCVal(2, top_mask) = delta_u;
FIndices(2, top_mask) = true;BCIndices = BCIndices(:);BCVal = BCVal(:);FIndices = FIndices(:);
%% 组装前准备
nQuad = 3; nQuadBdry = 2; max_iter = 1e4; 
Psi_plus_old = zeros(nQuad,size(IEN,2)); Psi_plus_rec = zeros(nQuad,size(IEN,2)); % 历史变量
node_range = 0:(nNodes-1);u_indices = (1:nNodes*nDoF)';u_indices(mod(u_indices-1, nDoF) >= nDim) = [];
d_indices = (node_range * nDoF) + (nDim + 1);BCIndices_u = BCIndices(u_indices);BCVal_u = BCVal(u_indices);
BCIndices_d = BCIndices(d_indices); BCVal_d = BCVal(d_indices); Sol_u = zeros(length(u_indices),1); Sol_d = zeros(length(d_indices),1);
active_indices_u = ~BCIndices_u; active_indices_d = ~BCIndices_d;
LM_u = reshape(repmat(reshape((IEN-1)*nDim, 1, []), nDim, 1)+ repmat((1:nDim)', 1, numel(IEN)), nNodesElement * nDim, []);  LM_d = IEN;
%% 加载
record_times_Sc = zeros(length(load_steps),6); Ru_tol = 1e-5; Rd_tol = 1e-5;
for step_no = 1:length(load_steps)
    run_time = 0;
    tic
    scu_times = 0;
    ntu_times = 0;
    scd_times = 0;
    ntd_times = 0;
    converged = false;
    Sol_u(BCIndices_u) = BCVal_u(BCIndices_u) * load_steps(step_no);
    Sol_d(BCIndices_d) = BCVal_d(BCIndices_d) * load_steps(step_no);
    for iter_no = 1:max_iter  %-------->Stepping for Newton iteration        
        if iter_no == 1
            [Ru, Psi_plus_rec] = Assemble_Ru(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
                constitutive, body_force, traction, load_steps(step_no), Psi_plus_rec, Psi_plus_old);
            Rd = Assemble_Rd(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
                constitutive, body_force, traction, load_steps(step_no), Psi_plus_rec, Psi_plus_old);
            Ru_active = Ru(active_indices_u);
            Rd_active = Rd(active_indices_d);
            Ru_ref = norm(Ru_active);
            Rd_ref = norm(Rd_active);
            record_times_Sc(step_no,1) =  load_steps(step_no);
        end
        
        [Sol_u,Sol_d,u_sctimes,u_nttimes,d_nttimes,~,doubleiter,~, Psi_plus_rec]=doublefieldschur_Neumann_space(Sol_u,Sol_d,Coord, IEN, LM_u, LM_d, elementType, ...
                    constitutive,body_force,traction, load_steps(step_no) ,step_no,'ud',Ru_ref,Rd_ref,...
                    active_indices_u,active_indices_d,Ru_tol,Rd_tol,Psi_plus_rec, Psi_plus_old);
                scu_times = scu_times + u_sctimes;
                ntu_times = ntu_times + u_nttimes;
                ntd_times = ntd_times + d_nttimes;
                converged = true;break
    end    
    if (converged)
        run_time = toc
        [Ru, ~] = Assemble_Ru(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
        constitutive, body_force, traction, load_steps(step_no), Psi_plus_rec, Psi_plus_old);
        force = sum(Ru(FIndices(:)==1));
        Psi_plus_old = Psi_plus_rec;
        record_times_Sc(step_no,3:9) = [scu_times, ntu_times, ntd_times, run_time, max(Sol_d), force, doubleiter-1];
        fprintf('Loading of %g calculation is complete.\n',load_steps(step_no));
%         fprintf('The displacement field and phase field have iterated %g and %g times, respectively.\n',u_times,d_times);
    else
        break;
    end
%% 保存信息
Sol(u_indices) = Sol_u;
Sol(d_indices) = Sol_d;
Answer=reshape(Sol, 3, []);
Displace=Answer(1:2,:);
D=Answer(3,:);
X = [Coord(1,IEN(1,:)); Coord(1,IEN(2,:)); Coord(1,IEN(3,:))];
Y = [Coord(2,IEN(1,:)); Coord(2,IEN(2,:)); Coord(2,IEN(3,:))];
D_elem = [Displace(1,IEN(1,:)); Displace(1,IEN(2,:)); Displace(1,IEN(3,:))];
figure(4)%
colormap('jet')
X = [Coord(1,IEN(1,:)); Coord(1,IEN(2,:)); Coord(1,IEN(3,:))];
Y = [Coord(2,IEN(1,:)); Coord(2,IEN(2,:)); Coord(2,IEN(3,:))];
D_elem = [D(IEN(1,:)); D(IEN(2,:)); D(IEN(3,:))];
XYZ=patch(X,Y,D_elem);
axis equal;
axis off;
shading interp;
colorbar;
caxis([0 1]);
saveas(XYZ, fullfile(outputfile, ['Phasefield-' strrep(num2str(load_steps(step_no), '%.5f'),'.','_') '.jpg']));
save(fullfile(outputfile, ['PF-' strrep(num2str(load_steps(step_no), '%.5f'),'.','_') '.mat']));
end

