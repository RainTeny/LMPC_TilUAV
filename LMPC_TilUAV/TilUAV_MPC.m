%% 针对倾转旋翼弹式轨迹跟踪问题的LMPC建模与控制
clear all
close all
t0=cputime;
%% 第一部分：Quadcopter state-space model
%% 系统参数
m=1; % mass kg
lm=0.18;%机臂长度0.18m
d=1.574e-7; % 单位N·m/(rad/s)^2，电机反扭距系数
g=9.81; 
Ix=7.5E-3; % 转动惯量 x 轴 kg*m^2
Iy=7.5E-3; % 转动惯量 y 轴 kg*m^2
Iz=1.3E-2; % 转动惯量 z 轴 kg*m^2
dx=0; dy=0; dz=0; % x y z 方向扰动为零
kx=0.25; ky=0.25; kz=0.25; % 空气阻力

%% 状态空间模型
%未简化的模型
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12
syms u1 u2 u3 u4 u5 u6
% x=[x y z vx vy vz roll pitch yaw wr wp wy]
% u=[sinα1*F1、sinα2*F3、cosα1*F1、cosα2*F3、F2、F4 ] F1-4分别是四个机翼升力
xdot1=x4;
xdot2=x5;
xdot3=x6;
xdot4=(u1+u2)*(cos(x8)*cos(x9))/m+(u3+u4+u5+u6)*(cos(x7)*sin(x8)*cos(x9)+sin(x7)*sin(x9))/m-kx*x4/m;
xdot5=(u1+u2)*(cos(x8)*sin(x9))/m+(u3+u4+u5+u6)*(cos(x7)*sin(x8)*cos(x9)-sin(x7)*cos(x9))/m-ky*x5/m;
xdot6=-(u1+u2)*sin(x8)/m+(u3+u4+u5+u6)*(cos(x7)*cos(x8))/m-g-kz*x6/m;
% xdot4=(u1/m)*(cos(x9)*sin(x8)*cos(x7)+sin(x9)*sin(x7))-kx*x4/m+dx/m;
% xdot5=(u1/m)*(sin(x9)*sin(x8)*cos(x7)-cos(x9)*sin(x7))-ky*x5/m+dy/m;
% xdot6=(u1/m)*cos(x8)*cos(x7)-g-kz*x6/m+dz/m;
xdot7=x10+x11*sin(x7)*tan(x8)+x12*cos(x7)*tan(x8);
xdot8=x11*cos(x7)-x12*sin(x7);
xdot9=sin(x7)*x11/cos(x8)+cos(x7)*x12/cos(x8);
xdot10=(lm*(u3-u4)/Ix)-((Iy-Iz)/Ix)*x11*x12;
xdot11=(lm*(u5-u6)/Iy)-((Iz-Ix)/Iy)*x10*x12;
xdot12=(lm*(u1+u2)/Iz)+(d*(u3+u4-u5-u6)/Iz)-((Ix-Iy)/Iz)*x10*x11;
% xdot10=(u2/Ix)-((Iy-Iz)/Ix)*x11*x12;
% xdot11=(u3/Iy)-((Iz-Ix)/Iy)*x10*x12;
% xdot12=(u4/Iz)-((Ix-Iy)/Iz)*x10*x11;
xdot=[xdot1,xdot2,xdot3,xdot4,xdot5,xdot6,xdot7,xdot8,xdot9,xdot10,xdot11,xdot12]';
x=[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12]';
u=[u1,u2,u3,u4,u5,u6]';
y=[x1,x2,x3,x7,x8,x9]';

%% 简化后的模型
% 为避免奇异化，使得雅各比矩阵求导不出来，近似看做 sin x -> x; 1/cos x -> 1; tan x -> x; no distrurbances
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12
syms u1 u2 u3 u4 u5 u6
% x=[x y z vx vy vz roll pitch yaw wr wp wy]
xdot1=x4;
xdot2=x5;
xdot3=x6;
xdot4=(u1+u2)*(cos(x8)*cos(x9))/m+(u3+u4+u5+u6)*(cos(x7)*sin(x8)*cos(x9)+sin(x7)*sin(x9))/m-kx*x4/m;
xdot5=(u1+u2)*(cos(x8)*sin(x9))/m+(u3+u4+u5+u6)*(cos(x7)*sin(x8)*cos(x9)-sin(x7)*cos(x9))/m-ky*x5/m;
xdot6=-(u1+u2)*sin(x8)/m+(u3+u4+u5+u6)*(cos(x7)*cos(x8))/m-g-kz*x6/m;
% xdot7=x10+x11*sin(x7)*tan(x8)+x12*cos(x7)*tan(x8);
xdot7=x10+x11*sin(x7)*x8+x12*cos(x7)*x8;
xdot8=x11*cos(x7)-x12*sin(x7);
% xdot9=sin(x7)*x11/cos(x8)+cos(x7)*x12/cos(x8);
xdot9=sin(x7)*x11/1+cos(x7)*x12/1
xdot10=(lm*(u3-u4)/Ix)-((Iy-Iz)/Ix)*x11*x12;
xdot11=(lm*(u5-u6)/Iy)-((Iz-Ix)/Iy)*x10*x12;
xdot12=(lm*(u1+u2)/Iz)+(d*(u3+u4-u5-u6)/Iz)-((Ix-Iy)/Iz)*x10*x11;
xdot=[xdot1 xdot2 xdot3 xdot4 xdot5 xdot6 xdot7 xdot8 xdot9 xdot10 xdot11 xdot12].';
x=[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12].';
u=[u1,u2,u3,u4,u5,u6].';
y=[x1,x2,x3,x9].';
[x_size,aux]=size(x);
[y_size,aux]=size(y);
[u_size,aux]=size(u);

%% 雅可比线性化
% 稳定点
ue=[0,0,0.25*m*g,0.25*m*g,0.25*m*g,0.25*m*g].';
xe=[x1,x2,x3,0,0,0,0,0,0,0,0,0].';%状态量初始值
JA=jacobian(xdot,x.');
JB=jacobian(xdot,u.');
JC=jacobian(y, x.');
A = double(subs(JA, [x; u], [xe; ue]));  % 带入稳定点，使用 double() 转数值】
B = double(subs(JB, [x; u], [xe; ue]));
C = double(subs(JC, [x; u], [xe; ue]));
%%离散状态空间模型
fs=50; % 50-100 Hz
Ts=1/fs; % 采样周期 = 0.02秒
sysc=ss(A,B,C,0);  %%创建实数或复数的状态空间模型/拉氏变换
sysd=c2d(sysc,Ts); %%将s域的表达式转化成z域的表达式/离散化 Zero-Order Hold（ZOH）
Am=sysd.A;
Bm=sysd.B;
Cm=sysd.C;
Dm=zeros(x_size,1);

%% MPC控制器
%% 初始化
N_sim=10*fs; %%模拟时间步数 samples = seconds*frequency
Nc=15;  % 控制时域（优化器只输出前 Nc 步的控制量）
Np=25; % 预测时域（输出预测跨度）
R = 0.001;   %控制输入增量的惩罚项权重（代价函数）
H = 1;
% Init control, reference and output signal
u=zeros(u_size,1);  
y=zeros(y_size,1); 

rx = zeros(1, N_sim);
ry = zeros(1, N_sim);
rz = zeros(1, N_sim);
ryaw = zeros(1, N_sim);
% ones(1,N_sim) ：此函数会生成一个长度为 N_sim 的行向量，其每个元素的值均为 1。
% 0*ones(1,N_sim) ：这是将全 1 向量与 0 相乘，其结果是得到一个长度为 N_sim 且元素全为 0 的行向量。
xm_vector=[];
u_vector=[];
deltau_vector=[];
y_vector=[];
%xm_vector、u_vector、deltau_vector 和 y_vector ：这些变量被初始化为空矩阵。

% 初始系统状态 全部为0
xm=zeros(x_size,1); % 状态向量
Xf=zeros(x_size+y_size,1); % 增强增量状态 [deltax y]'

%% 得到增广增量模型和增量控制参数
[GG, GF, GR, A, B, C , GD, F, G]=mpcgain_mimo(Am,Bm,Cm,Nc,Np,Dm); 

%% 增量轨迹控制的恒定部分
f1=H*GG+R*eye(Nc*u_size,Nc*u_size); % E for the cost function J=xEx'+x'F

%% 约束矩阵 M
ymax=[0.8 1 2];  % 输出最大值限制(1 x n_y)
ymin=[-0.3 -1 -1.5];  % 输出最小限制 (1 x n_y)
deltaumax=500;%控制量增量最大限制
deltaumin=-500;%控制量增量最小限制
umax=500;%控制量最大限制
umin=-500;%控制量最小限制

M_output=[];
[gm,gn]=size(G);
aux=eye(gn);
M_deltau=[aux(1,:);-aux(1,:);%提取第一行
    aux(2,:);-aux(2,:);
    aux(3,:);-aux(3,:);
    aux(4,:);-aux(4,:)];

M_u=[aux(1,:);-aux(1,:);
    aux(2,:);-aux(2,:);
    aux(3,:);-aux(3,:);
    aux(4,:);-aux(4,:)];

M=[M_output; M_deltau; M_u];%总约束矩阵
%alpha1_vector 和 alpha2_vector 是 global 变量或 workspace 中的变量，这样 TilUAV_plot 能读取
alpha1_vector = zeros(1, N_sim);
alpha2_vector = zeros(1, N_sim);

%% 循环
for k=1:N_sim 
    %% 参考轨迹
    rx(k)=sin(0.03*k);    
    ry(k)=cos(0.03*k);       
    rz(k)=0.1*k; 
%    rz(k)=sin(0.03*k);
    ryaw(k)=45*pi/180;
    
    r=[rx; ry; rz; ryaw];
    
    %% 计算控制输出 DeltaU
    %求第二部分 DeltaU    
    f2=GR*r(:,k)-GF*Xf; % F 表示成本函数 J=xEx'+x'F
    DeltaU=inv(f1)*f2; % 没有约束的控制输出 DeltaU 
    
    %% 没有约束下计算控制输出 
    gamma_output=[];
    gamma_deltau=[deltaumax;-deltaumin;
        deltaumax;-deltaumin;
        deltaumax;-deltaumin;
        deltaumax;-deltaumin];
    gamma_u=[umax-u(1);0+u(1);
        umax-u(2);-umin+u(2);
        umax-u(3);-umin+u(3);
        umax-u(4);-umin+u(4)];
    gamma=[gamma_output; gamma_deltau; gamma_u];
        
    %% 计算 u
    deltau=DeltaU(1:u_size);  % 应用滚动优化策略
    u=u+deltau; % (u(k)=u(k-1)+deltau(k))
    u = u - ue;  % 去除稳态输入（相对控制）
    %% 增量状态向量预测     
    xm_old=xm; 
    xm=Am*xm+Bm*u; 
    y=Cm*xm;    
    Xf=[xm-xm_old;y]; 
    
    %% 添加到向量
    xm_vector=[xm_vector xm];
    y_vector=[y_vector y];
    u_vector=[u_vector u];
    deltau_vector=[deltau_vector deltau];
  % 计算当前时间步的 alpha1 和 alpha2
    [a1, a2,U_1, U_2, U_3, U_4] = Til_inputU(u(1), u(2), u(3), u(4), u(5), u(6));
     alpha1_vector(k) = a1;
     alpha2_vector(k) = a2;

end
[alpha_1, alpha_2, U_1, U_2, U_3, U_4]=Til_inputU(u1,u2,u3,u4,u5,u6);
TilUAV_plot;

% 用于计算累计时间误差
errorsum=sum(sqrt((y_vector(1,:)-r(1,:)).^2)+sqrt((y_vector(2,:)-r(2,:)).^2)+sqrt((y_vector(3,:)-r(3,:)).^2)+sqrt((y_vector(4,:)-r(4,:)).^2))
% 计算时间
elapsed_time=cputime-t0