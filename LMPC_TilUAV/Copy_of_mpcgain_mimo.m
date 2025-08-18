function [Phi_Phi,Phi_F,Phi_R,A_e,B_e,C_e,Phi_D,F,Phi] = Copy_of_mpcgain_mimo(A,B,H,Nc,Np,D)
% mpcgain_mimo - 计算MPC控制的增益矩阵
% 输入:
%   A  - 系统状态矩阵
%   B  - 系统输入矩阵
%   H  - 输出矩阵C
%   Nc - 控制时域
%   Np - 预测时域
%   D  - 干扰矩阵
% 输出:
%   Phi_Phi - 二次规划中的H矩阵 (GG)
%   Phi_F   - 二次规划中的状态相关一次项 (GF)
%   Phi_R   - 二次规划中的参考相关一次项 (GR)
%   A_e     - 增广系统矩阵
%   B_e     - 增广输入矩阵
%   C_e     - 增广输出矩阵
%   Phi_D   - 干扰项系数矩阵
%   F       - 自由响应矩阵
%   Phi     - 系统响应矩阵(Block Toeplitz)

% 获取矩阵维度
[m1,n1] = size(H);    % H矩阵的行数m1和列数n1
[nb,n_in] = size(B);  % B矩阵的行数nb和列数n_in

% 检查D矩阵维度并修正
if size(D,1) ~= nb || size(D,2) ~= n_in
    D = zeros(nb, n_in);
    disp('警告: D矩阵维度已调整以匹配系统');
end
[nd,nd_in] = size(D);  % 修正后的D矩阵维度

% 创建增广状态空间矩阵
% 状态方程: x(k+1) = A*x(k) + B*u(k)
% 输出方程: y(k) = H*x(k)
% 增广系统: [x(k+1); y(k+1)] = [A, 0; H*A, I] * [x(k); y(k)] + [B; H*B] * u(k)

% 创建增广系统矩阵A_e = [A, 0; H*A, I]
A_e = zeros(nb+m1, nb+m1);
A_e(1:nb, 1:nb) = A;
A_e(nb+1:nb+m1, 1:nb) = H*A;
A_e(nb+1:nb+m1, nb+1:nb+m1) = eye(m1);

% 创建增广控制矩阵B_e = [B; H*B]
B_e = zeros(nb+m1, n_in);
B_e(1:nb, 1:n_in) = B;
B_e(nb+1:nb+m1, 1:n_in) = H*B;

% 创建增广干扰矩阵D_e = [D; H*D]
D_e = zeros(nb+m1, nd_in);
D_e(1:nb, 1:nd_in) = D;
D_e(nb+1:nb+m1, 1:nd_in) = H*D;

% 创建增广输出矩阵C_e = [0, I]，只关心增广状态的输出部分
C_e = zeros(m1, nb+m1);
C_e(1:m1, nb+1:nb+m1) = eye(m1);

% 计算自由响应矩阵F
F = [];
for i = 1:Np
    F = [F; C_e * (A_e^i)];
end

% 构建系统响应矩阵Phi (Block Toeplitz矩阵)
Phi = zeros(m1*Np, n_in*Nc);
for i = 1:Np
    for j = 1:min(i,Nc)
        Phi((i-1)*m1+1:i*m1, (j-1)*n_in+1:j*n_in) = C_e * (A_e^(i-j)) * B_e;
    end
end

% 计算干扰矩阵的影响
Phi_d = zeros(m1*Np, nd_in);
for i = 1:Np
    Phi_d((i-1)*m1+1:i*m1, :) = C_e * (A_e^(i-1)) * D_e;
end

% 计算二次规划的系数矩阵
Phi_Phi = Phi' * Phi;        % 二次项系数
Phi_F = Phi' * F;            % 一次项系数(状态相关)
Phi_R = Phi' * eye(m1*Np);   % 一次项系数(参考相关)

% 调整Phi_R维度以匹配参考向量
Phi_R = zeros(n_in*Nc, m1);
for i = 1:Np
    Phi_R = Phi_R + Phi(:, 1:n_in*Nc)' * eye(m1*Np, m1*Np) * eye(m1*Np, m1);
end

% 计算干扰项系数
if nd_in > 0
    Phi_D = Phi' * Phi_d;
else
    Phi_D = zeros(n_in*Nc, 1);
end
end