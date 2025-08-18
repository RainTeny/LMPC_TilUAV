%% 此函数用于得到增广增量模型和增量控制参数
function [Phi_Phi,Phi_F,Phi_R,A_e,B_e,C_e,Phi_D,F,Phi] = mpcgain_mimo(A,B,H,Nc,Np,D)
% Nc控制时域,Np预测时域，H=输出矩阵C
global data

[m1,n1]=size(H);% 获取矩阵 H 的行数 m1 和列数 n1
[nb,n_in]=size(B); % 获取矩阵 B 的行数 nb 和列数 n_in
[nd,nd_in]=size(D);%未使用
[y1 y2]=size(H*B);
[q,q1]=size(A);
A_e=zeros(q+m1,q+m1); % 创建一个大小为 (q+m1)×(q+m1) 的零矩阵 A_e
A_e(1:q,1:q)=A;
A_e(q+1:q+m1,1:n1)=H*A; % 将矩阵乘积 H*A 放置在 A_e 的左下角区域
A_e(q+1:m1+q,q+1:q+m1)=eye(m1); % 将 m1×m1 的单位矩阵放置在 A_e 的右下角区域
% A_e = [A, 0;
%        C*A, eye(p)];
% B_e = [B;
%        C*B];
% C_e = [0, I];


% B_e(1:nb,1:n_in)=B
B_e=B;
B_e(nb+1:nb+y1,1:y2)=H*B;

[y1 y2]=size(H*D);
D_e=D;
D_e(nd+1:nd+y1,1:y2)=H*D;

C_e=zeros(m1,n1);
C_e(1:m1,n1+1:n1+m1)=eye(m1);


[x1,x2]=size(C_e*A_e);
for kk=1:Np
    nn=kk-1;
F(x1*nn+1:x1*nn+x1,1:x2)=C_e*A_e^kk;%系统自由响应矩阵F
end

[x3,x4]=size(C_e*B_e);% 输出维度 × 控制输入维度
for i=1:Nc
    mm=i-1;
    for j=1:Np
        nn=j-1;
        if j<i
         Phi(x3*nn+1:x3*nn+x3,x4*mm+1:x4*mm+x4)=zeros(x3,x4);
        else
          Phi(x3*nn+1:x3*nn+x3,x4*mm+1:x4*mm+x4)=C_e*A_e^(j-i)*B_e;
        end
    end
end%控制输入对输出的影响矩阵 Φ， 下三角的 Block Toeplitz 矩阵

[x3,x4]=size(C_e*D_e);
for i=1:Nc
    mm=i-1;
    for j=1:Np
        nn=j-1;
        if j<i
         Phi_d(x3*nn+1:x3*nn+x3,x4*mm+1:x4*mm+x4)=zeros(x3,x4);
        else
          Phi_d(x3*nn+1:x3*nn+x3,x4*mm+1:x4*mm+x4)=C_e*A_e^(j-i)*D_e;
        end
    end
end
% 干扰D对输出的影响矩阵 Φ
[n,m]=size(C_e);
% W=[eye(m) zeros(m,Nc)];
% BarRs=eye(n,m+n);%[0 1 0 1 0 1 0 1]';%ones(Np*2,1);%1
[x1,x2]=size(Phi);
Phi_Phi = Phi' * Phi;  % 二次项（ΔU^T * ΦᵀΦ * ΔU）

% 在 Phi_Phi 中增加位置误差的惩罚项
% Phi_Phi = Phi_Phi + 0.01 * eye(size(Phi_Phi));  % 通过加大对位置的惩罚来增加对误差的惩罚

Phi_F = Phi' * F;      % 一次项中的 Xf 项（ΔU^T * ΦᵀF * Xf）
Phi_D = Phi' * Phi_d;  % 干扰项（ΔU^T * Φᵀ * Φd * d）
[x1,x2]=size(Phi_F);
Phi_R=Phi_F(1:x1,x2-m1+1:x2);
[x1,x2]=size(Phi_D);
Phi_D=10*Phi_D(:,1);
% Phi_D=Phi_D(1:x1,x2-m1+1:x2);
end