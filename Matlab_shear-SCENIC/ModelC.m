function [StressVector, D, Psi_plus, sigma_plus, gc, ell, gd_diff] = ModelC(StrainVector, d)
E = 210e3;
nu = .3;
gc = 2.7;
sigma_c = 2445.42;
ell = 0.01;
a1=27/128*E*gc/ell/sigma_c^2;
%%
lambda = E*nu/((1+nu)*(1-2*nu));
mu = E/(2*(1+nu));
kmin = 1e-10 * E;
%% 
gd = (1-d)^2+1e-8;
gd_diff = -2*(1-d);
% %% Wu
% gd = (1-d)^2/((1-d)^2+a1*d*(1-0.5*d))+1e-8;
% gd_diff = -a1*(1-d)/((1-d)^2+a1*d*(1-0.5*d))^2;

angle_bracket_Plus = @(a) (a+abs(a))/2;
angle_bracket_Minus = @(a) (a-abs(a))/2;
switch size(StrainVector, 1)
    case 3
        nDim = 2;
        StrainTensor = [StrainVector(1), StrainVector(3)/2, 0
            StrainVector(3)/2, StrainVector(2), 0
            0, 0, 0];
    case 6
        nDim = 3;
        StrainTensor = [StrainVector(1), StrainVector(6)/2, StrainVector(5)/2
            StrainVector(6)/2, StrainVector(2), StrainVector(4)/2
            StrainVector(5)/2, StrainVector(4)/2, StrainVector(3)];
end

[principal_directions, principal_strains] = eig(StrainTensor);
%注意，principal_directions的第�?列才是向量n1，n12代表向量a1的第二个分量
%和组内公众号相互呼应
n11 = principal_directions(1,1);
n12 = principal_directions(2,1);
n13 = principal_directions(3,1);
n21 = principal_directions(1,2);
n22 = principal_directions(2,2);
n23 = principal_directions(3,2);
n31 = principal_directions(1,3);
n32 = principal_directions(2,3);
n33 = principal_directions(3,3);
%A_eps代表应变的转置矩阵，可以参�?�公众号文章https://mp.weixin.qq.com/s/rORT3N3B4z2PK-4bEQyASQ
%或�?�组内implementation details for miehe....pdf
A_eps = [n11^2 n12^2 n13^2 n12*n13 n11*n13 n11*n12; 
         n21^2 n22^2 n23^2 n22*n23 n21*n23 n21*n22; 
         n31^2 n32^2 n33^2 n32*n33 n31*n33 n31*n32; 
         2*n21*n31 2*n22*n32 2*n23*n33 n22*n33+n23*n32 n21*n33+n23*n31 n21*n32+n22*n31; 
         2*n11*n31 2*n12*n32 2*n13*n33 n12*n33+n13*n32 n11*n33+n13*n31 n11*n32+n12*n31; 
         2*n11*n21 2*n12*n22 2*n13*n23 n12*n23+n13*n22 n11*n23+n13*n21 n11*n22+n12*n21];
principal_strains = diag(principal_strains);
%4.eig函数求出来的特征方向，特征�?�均是以矩阵形式输出的？ ANS:yes，所以我们接下来
%在表示特征�?�的时�?�，继续用了diag函数把对角线元素取出来�??
Psi_plus = lambda/2 * angle_bracket_Plus(trace(StrainTensor))^2 + ...
    mu * sum(angle_bracket_Plus(principal_strains).^2);  %此处利用sum函数求和非常�?
    lambdaHeavplus = lambda*heaviside(trace(StrainTensor));
    lambdaHeavminus = lambda*heaviside(-trace(StrainTensor));
Dn_Plus = [lambdaHeavplus+2*mu*heaviside(principal_strains(1)) ...
    lambdaHeavplus lambdaHeavplus 0 0 0; 
    lambdaHeavplus lambdaHeavplus+2*mu*heaviside(principal_strains(2)) ...
    lambdaHeavplus 0 0 0; 
    lambdaHeavplus lambdaHeavplus ...
    lambdaHeavplus+2*mu*heaviside(principal_strains(3)) 0 0 0 ; 
    0 0 0 mu*Hab_Plus(principal_strains(2),principal_strains(3)) 0 0; 
    0 0 0 0 mu*Hab_Plus(principal_strains(1),principal_strains(3)) 0; 
    0 0 0 0 0 mu*Hab_Plus(principal_strains(1),principal_strains(2))];
Dn_Minus = [lambdaHeavminus+2*mu*heaviside(-principal_strains(1)) ...
    lambdaHeavminus lambdaHeavminus 0 0 0; ...
    lambdaHeavminus lambdaHeavminus+2*mu*heaviside(-principal_strains(2)) ...
    lambdaHeavminus 0 0 0; ...
    lambdaHeavminus lambdaHeavminus ...
    lambdaHeavminus+2*mu*heaviside(-principal_strains(3)) 0 0 0 ; ...
    0 0 0 mu*Hab_Minus(principal_strains(2),principal_strains(3)) 0 0; ...
    0 0 0 0 mu*Hab_Minus(principal_strains(1),principal_strains(3)) 0; ...
    0 0 0 0 0 mu*Hab_Minus(principal_strains(1),principal_strains(2))];
D_Plus = (A_eps)' * Dn_Plus * A_eps;
D      = A_eps' * (gd*Dn_Plus + Dn_Minus) * A_eps;
switch nDim
    case 3        
        %转化成张�?/矩阵的形式时，需注意：对�?2D情况，StressVector�?3）代表\sigma_{12}�?
        %与应变向量和张量之间的转换关系不�?样，应变向量的StrainVector�?3）代�?2*epsilon_{12}
        Stress_plus = D_Plus * StrainVector;
        sigma_plus = [Stress_plus(1) Stress_plus(6) Stress_plus(5); ...
            Stress_plus(6) Stress_plus(2) Stress_plus(4); ...
            Stress_plus(5) Stress_plus(4) Stress_plus(3)];
    case 2
        D = D([1,2,6], [1,2,6]);
        D_Plus = D_Plus([1,2,6], [1,2,6]);
        StressVector = D * StrainVector;
        sigma_plus = D_Plus * StrainVector;
end
end
function [s] =Hab_Plus(a,b)
angle_bracket_Plus = @(a) (a+abs(a))/2;
if a==b 
    s=heaviside(a);
else
    s=(angle_bracket_Plus(a) - angle_bracket_Plus(b))/(a-b);
end
end
function [s] =Hab_Minus(a,b)
angle_bracket_Minus = @(a) (a-abs(a))/2;
if a==b 
    s=heaviside(-a);
else
    s=(angle_bracket_Minus(a) - angle_bracket_Minus(b))/(a-b);
end
end



