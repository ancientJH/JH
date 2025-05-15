function [Sol_u,Sol_d,a_sctimes,a_nttimes,b_nttimes, rec_Pd, doubleiter, rec_R, Psi_plus_rec] = doublefieldschur_Neumann_space(Sol_u,Sol_d,Coord, IEN, LM_u, LM_d, elementType, ...
    constitutive,body_force,traction, load_pre ,step_no,whichhalf,Ru_ref,Rd_ref,...
    active_indices_u,active_indices_d,Ru_tol,Rd_tol,Psi_plus_rec, Psi_plus_old)
flag = 's求解';
open_eta = 0;%是否采用eta
rec_Rus = [0,0,0];
rec_Rusc = [0,0,0];
rec_Pd = [0,0,0,0];
rho_ref = 0.25;
rec_R = [];

a_nttimes = 0;
a_sctimes = 0;
b_nttimes = 0;
% eta_max = 1;
% delta_1=zeros(sum(active_indices_u),1);
M=cell(3,1);
for doubleiter=1:500

    [Ru, Psi_plus_rec] = Assemble_Ru(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
        constitutive, body_force, traction, load_pre, Psi_plus_rec, Psi_plus_old);
    Rd = Assemble_Rd(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
        constitutive, body_force, traction, load_pre, Psi_plus_rec, Psi_plus_old);
    Ru_active = Ru(active_indices_u);
    Rd_active = Rd(active_indices_d);
    %     rec_norm = [rec_norm;norm(Ru_active)];
    if strcmp(flag,'s求解')
        rec_Rus = [rec_Rus,norm(Ru_active)];rec_Rus = rec_Rus(:,2:4);
        if (rec_Rus(3)/rec_Rus(2)>rho_ref) && (rec_Rus(3)*rec_Rus(2)~=0)
            flag='sc法';
        end
    else
        rec_Rusc = [rec_Rusc,norm(Ru_active)];rec_Rusc = rec_Rusc(:,2:4);
        if ((rec_Rusc(3)/rec_Rusc(2)>1)||(rec_Rusc(2)/rec_Rusc(1)>1))&&(rec_Rusc(1)*rec_Rusc(2)*rec_Rusc(3)~=0)
            open_eta = 1;
        end
    end

    fprintf('load %g, iter. %d, L2-norm of Ru = %g, L2-norm of Rd = %g\n',load_pre, doubleiter-1,norm(Ru_active),norm(Rd_active));

    rec_R = [rec_R;norm(Ru_active)];
    if norm(Ru_active)<= Ru_tol * Ru_ref&&norm(Rd_active)<= Rd_tol * Rd_ref
        break
    end

    if strcmp(flag,'s求解')
        k=0;
        %牛顿法迭代1至收敛
        nt1_converged = false;
        for nt_1=1:100
            Ku = Assemble_Ku(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
                constitutive, body_force, traction, load_pre, Psi_plus_rec, Psi_plus_old);
            Kuu_active = Ku(active_indices_u, active_indices_u);
            Sol_u(active_indices_u) = Sol_u(active_indices_u) - Kuu_active\Ru_active;

            a_nttimes = a_nttimes + 1;
            [Ru, Psi_plus_rec] = Assemble_Ru(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
                constitutive, body_force, traction, load_pre, Psi_plus_rec, Psi_plus_old);
            Ru_active = Ru(active_indices_u);
            fprintf('L2-norm of Ru = %g\n',norm(Ru_active));
            if norm(Ru_active)<= Ru_tol * Ru_ref
                nt1_converged = true;break;
            end
        end
        if ~nt1_converged
            converged = false;
            error('牛顿法失败')
        end
    else

        [Ru, Psi_plus_rec] = Assemble_Ru(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
            constitutive, body_force, traction, load_pre, Psi_plus_rec, Psi_plus_old);
        Rd = Assemble_Rd(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
            constitutive, body_force, traction, load_pre, Psi_plus_rec, Psi_plus_old);
        [Kuu, Kdd, Kud, Kdu] = Assemble_K(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
            constitutive, body_force, traction, load_pre, Psi_plus_rec, Psi_plus_old);
        Ru_active = Ru(active_indices_u);
        Rd_active = Rd(active_indices_d);%这里需要改
        Kuu_active = Kuu(active_indices_u, active_indices_u);
        Kdd_active = Kdd(active_indices_d, active_indices_d);
        Kud_active = Kud(active_indices_u, active_indices_d);
        Kdu_active = Kdu(active_indices_d, active_indices_u);
        
        clear Kuu Kdd Kud Kdu
     
        if open_eta == 1
            zeros_uu = sparse(size(Kuu_active,1),size(Kuu_active,2));
            zeros_dd = sparse(size(Kdd_active,1),size(Kdd_active,2));
            zeros_ud = sparse(size(Kud_active,1),size(Kud_active,2));
            zeros_du = sparse(size(Kdu_active,1),size(Kdu_active,2));
            lambda = abs(eigs(-sparse([zeros_uu, Kud_active; zeros_du, zeros_dd]), sparse([Kuu_active,zeros_ud; Kdu_active, Kdd_active]), 1));
            if lambda<1
                eta = 1;
            else
                eta = sqrt(1/lambda)*0.9;
            end
            if lambda<0.25
                flag = 's求解';
            end
            if eta<0
                eta = 0;
            end


            delta_1 = approxd1(eta*Kud_active,Kdd_active,eta*Kdu_active,Kuu_active,Ru_active,Rd_active,50,1e-6);
            
            Sol_u(active_indices_u) = Sol_u(active_indices_u) + delta_1;
            %                     rec_eta = [rec_eta;eta];
            a_sctimes = a_sctimes + 1;

            fprintf('eta = %g\n',eta);
        else

            delta_1 = approxd1(Kud_active,Kdd_active,Kdu_active,Kuu_active,Ru_active,Rd_active,50,1e-6);


            Sol_u(active_indices_u) = Sol_u(active_indices_u) + delta_1;
            
            a_sctimes = a_sctimes + 1;
        end
        [~, Psi_plus_rec] = Assemble_Ru(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
            constitutive, body_force, traction, load_pre, Psi_plus_rec, Psi_plus_old);
    end
    Rd = Assemble_Rd(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
            constitutive, body_force, traction, load_pre, Psi_plus_rec, Psi_plus_old);
    Rd_active = Rd(active_indices_d);
    fprintf('L2-norm of Rd = %g\n',norm(Rd_active));
    if norm(Ru_active) <= Ru_tol * Ru_ref && norm(Rd_active) <= Rd_tol * Rd_ref
        converged = true;
        break;
    end
    nt2_converged = false;
    for nt_2=1:100
        Kdd = Assemble_Kd(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
                constitutive, body_force, traction, load_pre, Psi_plus_rec, Psi_plus_old);
        Kdd_active = Kdd(active_indices_d, active_indices_d);
        Sol_d(active_indices_d) = Sol_d(active_indices_d) - Kdd_active\Rd_active;

        b_nttimes = b_nttimes + 1;
        Rd = Assemble_Rd(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
            constitutive, body_force, traction, load_pre, Psi_plus_rec, Psi_plus_old);
        Rd_active = Rd(active_indices_d);
        fprintf('L2-norm of Rd = %g\n',norm(Rd_active));
        if norm(Rd_active)<= Rd_tol * Rd_ref
            nt2_converged = true;break;
        end
    end
    if ~nt2_converged
        converged = false;
        error('牛顿法失败')
    end

    if doubleiter == 500
        error('不收敛')
    end
end
end


function [delta_1, resvec] = approxd1(Kud_active, Kd_active, Kdu_active, Ku_active, Ru_active, Rd_active, max_iter, tol)
% 输入参数赋值
Kuu = Ku_active;
Kud = Kud_active;
Kdd = Kd_active;
Kdu = Kdu_active;
Ru = Ru_active;
Rd = Rd_active;

% 计算 y_i 和定义 A, b
y_i = Kuu \ (Ru - Kud * (Kdd \ Rd));
M = @(x) Kuu \ (Kud * (Kdd \ (Kdu * x)));
A = @(x) x - M(x);
b = -y_i;
n = length(b);

% 初始化
if nargin < 8, tol = 1e-6; end  % 默认容差
if nargin < 7, max_iter = 100; end  % 默认最大迭代次数
u0 = zeros(n, 1);  % 初始解
k = 1;
r0 = b - A(u0);  % 初始残差
res = norm(r0);  % 初始残差范数
rhs = zeros(max_iter + 1, 1);  % 初始化 rhs 为列向量
rhs(1) = res;  % 初始右端向量
Q = zeros(n, max_iter + 1);  % 预分配 Q
Q(:, 1) = r0 / res;  % 第一个正交向量
H = zeros(max_iter + 1, max_iter);  % Hessenberg 矩阵
R = zeros(max_iter + 1, max_iter);  % 上三角矩阵
c = zeros(max_iter, 1);  % Givens 旋转的余弦
s = zeros(max_iter, 1);  % Givens 旋转的正弦
resvec = zeros(max_iter + 1, 1);  % 残差历史
resvec(1) = res;

% GMRES 主循环
while resvec(k) / norm(b) > tol && k <= max_iter
    % Arnoldi 过程
    Q(:, k+1) = A(Q(:, k));  % 计算新向量
    for j = 1:k  % Gram-Schmidt 正交化
        H(j, k) = Q(:, k+1)' * Q(:, j);
        Q(:, k+1) = Q(:, k+1) - H(j, k) * Q(:, j);
    end
    H(k+1, k) = norm(Q(:, k+1));  % 计算 H 的下对角元素
    if H(k+1, k) ~= 0  % 避免除以零
        Q(:, k+1) = Q(:, k+1) / H(k+1, k);
    end

    % QR 分解（Givens 旋转）
    R(:, k) = H(:, k);  % 将 H 的第 k 列复制到 R
    for j = 1:k-1  % 应用之前的 Givens 旋转
        Rk = c(j) * R(j, k) + s(j) * R(j+1, k);
        R(j+1, k) = -s(j) * R(j, k) + c(j) * R(j+1, k);
        R(j, k) = Rk;
    end
    % 计算新的 Givens 旋转参数
    rho = sqrt(R(k, k)^2 + R(k+1, k)^2);
    if rho ~= 0  % 避免除以零
        c(k) = R(k, k) / rho;
        s(k) = R(k+1, k) / rho;
        R(k, k) = rho;
        R(k+1, k) = 0;
    end

    rhs(k+1) = -s(k) * rhs(k);
    rhs(k) = c(k) * rhs(k);
    resvec(k+1) = abs(rhs(k+1));  % 更新残差

    k = k + 1;
end
fprintf('第%d次迭代收敛，相对残差为%d\n',k,resvec(k) / norm(b))
if k >= max_iter
    warn('超出最大迭代次数')
end
% 求解
k = k - 1;  % 调整 k 到最后一轮迭代
y = R(1:k, 1:k) \ rhs(1:k);  % 求解上三角系统
delta_1 = u0 + Q(:, 1:k) * y;  % 计算最终解
end
