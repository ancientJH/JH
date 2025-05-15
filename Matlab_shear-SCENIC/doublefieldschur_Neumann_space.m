function [Sol_u,Sol_d,a_sctimes,a_nttimes,b_nttimes, rec_Pd, doubleiter, rec_R, Psi_plus_rec] = doublefieldschur_Neumann_space(Sol_u,Sol_d,Coord, IEN, LM_u, LM_d, elementType, ...
    constitutive,body_force,traction, load_pre ,step_no,whichhalf,Ru_ref,Rd_ref,...
    active_indices_u,active_indices_d,Ru_tol,Rd_tol,Psi_plus_rec, Psi_plus_old)
flag = 's���';
open_eta = 0;%�Ƿ����eta
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
    if strcmp(flag,'s���')
        rec_Rus = [rec_Rus,norm(Ru_active)];rec_Rus = rec_Rus(:,2:4);
        if (rec_Rus(3)/rec_Rus(2)>rho_ref) && (rec_Rus(3)*rec_Rus(2)~=0)
            flag='sc��';
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

    if strcmp(flag,'s���')
        k=0;
        %ţ�ٷ�����1������
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
            error('ţ�ٷ�ʧ��')
        end
    else

        [Ru, Psi_plus_rec] = Assemble_Ru(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
            constitutive, body_force, traction, load_pre, Psi_plus_rec, Psi_plus_old);
        Rd = Assemble_Rd(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
            constitutive, body_force, traction, load_pre, Psi_plus_rec, Psi_plus_old);
        [Kuu, Kdd, Kud, Kdu] = Assemble_K(Sol_u, Sol_d, Coord, IEN, LM_u, LM_d, elementType, ...
            constitutive, body_force, traction, load_pre, Psi_plus_rec, Psi_plus_old);
        Ru_active = Ru(active_indices_u);
        Rd_active = Rd(active_indices_d);%������Ҫ��
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
                flag = 's���';
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
        error('ţ�ٷ�ʧ��')
    end

    if doubleiter == 500
        error('������')
    end
end
end


function [delta_1, resvec] = approxd1(Kud_active, Kd_active, Kdu_active, Ku_active, Ru_active, Rd_active, max_iter, tol)
% ���������ֵ
Kuu = Ku_active;
Kud = Kud_active;
Kdd = Kd_active;
Kdu = Kdu_active;
Ru = Ru_active;
Rd = Rd_active;

% ���� y_i �Ͷ��� A, b
y_i = Kuu \ (Ru - Kud * (Kdd \ Rd));
M = @(x) Kuu \ (Kud * (Kdd \ (Kdu * x)));
A = @(x) x - M(x);
b = -y_i;
n = length(b);

% ��ʼ��
if nargin < 8, tol = 1e-6; end  % Ĭ���ݲ�
if nargin < 7, max_iter = 100; end  % Ĭ������������
u0 = zeros(n, 1);  % ��ʼ��
k = 1;
r0 = b - A(u0);  % ��ʼ�в�
res = norm(r0);  % ��ʼ�в��
rhs = zeros(max_iter + 1, 1);  % ��ʼ�� rhs Ϊ������
rhs(1) = res;  % ��ʼ�Ҷ�����
Q = zeros(n, max_iter + 1);  % Ԥ���� Q
Q(:, 1) = r0 / res;  % ��һ����������
H = zeros(max_iter + 1, max_iter);  % Hessenberg ����
R = zeros(max_iter + 1, max_iter);  % �����Ǿ���
c = zeros(max_iter, 1);  % Givens ��ת������
s = zeros(max_iter, 1);  % Givens ��ת������
resvec = zeros(max_iter + 1, 1);  % �в���ʷ
resvec(1) = res;

% GMRES ��ѭ��
while resvec(k) / norm(b) > tol && k <= max_iter
    % Arnoldi ����
    Q(:, k+1) = A(Q(:, k));  % ����������
    for j = 1:k  % Gram-Schmidt ������
        H(j, k) = Q(:, k+1)' * Q(:, j);
        Q(:, k+1) = Q(:, k+1) - H(j, k) * Q(:, j);
    end
    H(k+1, k) = norm(Q(:, k+1));  % ���� H ���¶Խ�Ԫ��
    if H(k+1, k) ~= 0  % ���������
        Q(:, k+1) = Q(:, k+1) / H(k+1, k);
    end

    % QR �ֽ⣨Givens ��ת��
    R(:, k) = H(:, k);  % �� H �ĵ� k �и��Ƶ� R
    for j = 1:k-1  % Ӧ��֮ǰ�� Givens ��ת
        Rk = c(j) * R(j, k) + s(j) * R(j+1, k);
        R(j+1, k) = -s(j) * R(j, k) + c(j) * R(j+1, k);
        R(j, k) = Rk;
    end
    % �����µ� Givens ��ת����
    rho = sqrt(R(k, k)^2 + R(k+1, k)^2);
    if rho ~= 0  % ���������
        c(k) = R(k, k) / rho;
        s(k) = R(k+1, k) / rho;
        R(k, k) = rho;
        R(k+1, k) = 0;
    end

    rhs(k+1) = -s(k) * rhs(k);
    rhs(k) = c(k) * rhs(k);
    resvec(k+1) = abs(rhs(k+1));  % ���²в�

    k = k + 1;
end
fprintf('��%d�ε�����������Բв�Ϊ%d\n',k,resvec(k) / norm(b))
if k >= max_iter
    warn('��������������')
end
% ���
k = k - 1;  % ���� k �����һ�ֵ���
y = R(1:k, 1:k) \ rhs(1:k);  % ���������ϵͳ
delta_1 = u0 + Q(:, 1:k) * y;  % �������ս�
end
