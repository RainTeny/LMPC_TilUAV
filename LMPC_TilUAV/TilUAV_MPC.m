%% �����ת����ʽ�켣���������LMPC��ģ�����
clear all
close all
t0=cputime;
%% ��һ���֣�Quadcopter state-space model
%% ϵͳ����
m=1; % mass kg
lm=0.18;%���۳���0.18m
d=1.574e-7; % ��λN��m/(rad/s)^2�������Ť��ϵ��
g=9.81; 
Ix=7.5E-3; % ת������ x �� kg*m^2
Iy=7.5E-3; % ת������ y �� kg*m^2
Iz=1.3E-2; % ת������ z �� kg*m^2
dx=0; dy=0; dz=0; % x y z �����Ŷ�Ϊ��
kx=0.25; ky=0.25; kz=0.25; % ��������

%% ״̬�ռ�ģ��
%δ�򻯵�ģ��
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12
syms u1 u2 u3 u4 u5 u6
% x=[x y z vx vy vz roll pitch yaw wr wp wy]
% u=[sin��1*F1��sin��2*F3��cos��1*F1��cos��2*F3��F2��F4 ] F1-4�ֱ����ĸ���������
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

%% �򻯺��ģ��
% Ϊ�������컯��ʹ���Ÿ��Ⱦ����󵼲����������ƿ��� sin x -> x; 1/cos x -> 1; tan x -> x; no distrurbances
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

%% �ſɱ����Ի�
% �ȶ���
ue=[0,0,0.25*m*g,0.25*m*g,0.25*m*g,0.25*m*g].';
xe=[x1,x2,x3,0,0,0,0,0,0,0,0,0].';%״̬����ʼֵ
JA=jacobian(xdot,x.');
JB=jacobian(xdot,u.');
JC=jacobian(y, x.');
A = double(subs(JA, [x; u], [xe; ue]));  % �����ȶ��㣬ʹ�� double() ת��ֵ��
B = double(subs(JB, [x; u], [xe; ue]));
C = double(subs(JC, [x; u], [xe; ue]));
%%��ɢ״̬�ռ�ģ��
fs=50; % 50-100 Hz
Ts=1/fs; % �������� = 0.02��
sysc=ss(A,B,C,0);  %%����ʵ��������״̬�ռ�ģ��/���ϱ任
sysd=c2d(sysc,Ts); %%��s��ı��ʽת����z��ı��ʽ/��ɢ�� Zero-Order Hold��ZOH��
Am=sysd.A;
Bm=sysd.B;
Cm=sysd.C;
Dm=zeros(x_size,1);

%% MPC������
%% ��ʼ��
N_sim=10*fs; %%ģ��ʱ�䲽�� samples = seconds*frequency
Nc=15;  % ����ʱ���Ż���ֻ���ǰ Nc ���Ŀ�������
Np=25; % Ԥ��ʱ�����Ԥ���ȣ�
R = 0.001;   %�������������ĳͷ���Ȩ�أ����ۺ�����
H = 1;
% Init control, reference and output signal
u=zeros(u_size,1);  
y=zeros(y_size,1); 

rx = zeros(1, N_sim);
ry = zeros(1, N_sim);
rz = zeros(1, N_sim);
ryaw = zeros(1, N_sim);
% ones(1,N_sim) ���˺���������һ������Ϊ N_sim ������������ÿ��Ԫ�ص�ֵ��Ϊ 1��
% 0*ones(1,N_sim) �����ǽ�ȫ 1 ������ 0 ��ˣ������ǵõ�һ������Ϊ N_sim ��Ԫ��ȫΪ 0 ����������
xm_vector=[];
u_vector=[];
deltau_vector=[];
y_vector=[];
%xm_vector��u_vector��deltau_vector �� y_vector ����Щ��������ʼ��Ϊ�վ���

% ��ʼϵͳ״̬ ȫ��Ϊ0
xm=zeros(x_size,1); % ״̬����
Xf=zeros(x_size+y_size,1); % ��ǿ����״̬ [deltax y]'

%% �õ���������ģ�ͺ��������Ʋ���
[GG, GF, GR, A, B, C , GD, F, G]=mpcgain_mimo(Am,Bm,Cm,Nc,Np,Dm); 

%% �����켣���Ƶĺ㶨����
f1=H*GG+R*eye(Nc*u_size,Nc*u_size); % E for the cost function J=xEx'+x'F

%% Լ������ M
ymax=[0.8 1 2];  % ������ֵ����(1 x n_y)
ymin=[-0.3 -1 -1.5];  % �����С���� (1 x n_y)
deltaumax=500;%�����������������
deltaumin=-500;%������������С����
umax=500;%�������������
umin=-500;%��������С����

M_output=[];
[gm,gn]=size(G);
aux=eye(gn);
M_deltau=[aux(1,:);-aux(1,:);%��ȡ��һ��
    aux(2,:);-aux(2,:);
    aux(3,:);-aux(3,:);
    aux(4,:);-aux(4,:)];

M_u=[aux(1,:);-aux(1,:);
    aux(2,:);-aux(2,:);
    aux(3,:);-aux(3,:);
    aux(4,:);-aux(4,:)];

M=[M_output; M_deltau; M_u];%��Լ������
%alpha1_vector �� alpha2_vector �� global ������ workspace �еı��������� TilUAV_plot �ܶ�ȡ
alpha1_vector = zeros(1, N_sim);
alpha2_vector = zeros(1, N_sim);

%% ѭ��
for k=1:N_sim 
    %% �ο��켣
    rx(k)=sin(0.03*k);    
    ry(k)=cos(0.03*k);       
    rz(k)=0.1*k; 
%    rz(k)=sin(0.03*k);
    ryaw(k)=45*pi/180;
    
    r=[rx; ry; rz; ryaw];
    
    %% ���������� DeltaU
    %��ڶ����� DeltaU    
    f2=GR*r(:,k)-GF*Xf; % F ��ʾ�ɱ����� J=xEx'+x'F
    DeltaU=inv(f1)*f2; % û��Լ���Ŀ������ DeltaU 
    
    %% û��Լ���¼��������� 
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
        
    %% ���� u
    deltau=DeltaU(1:u_size);  % Ӧ�ù����Ż�����
    u=u+deltau; % (u(k)=u(k-1)+deltau(k))
    u = u - ue;  % ȥ����̬���루��Կ��ƣ�
    %% ����״̬����Ԥ��     
    xm_old=xm; 
    xm=Am*xm+Bm*u; 
    y=Cm*xm;    
    Xf=[xm-xm_old;y]; 
    
    %% ��ӵ�����
    xm_vector=[xm_vector xm];
    y_vector=[y_vector y];
    u_vector=[u_vector u];
    deltau_vector=[deltau_vector deltau];
  % ���㵱ǰʱ�䲽�� alpha1 �� alpha2
    [a1, a2,U_1, U_2, U_3, U_4] = Til_inputU(u(1), u(2), u(3), u(4), u(5), u(6));
     alpha1_vector(k) = a1;
     alpha2_vector(k) = a2;

end
[alpha_1, alpha_2, U_1, U_2, U_3, U_4]=Til_inputU(u1,u2,u3,u4,u5,u6);
TilUAV_plot;

% ���ڼ����ۼ�ʱ�����
errorsum=sum(sqrt((y_vector(1,:)-r(1,:)).^2)+sqrt((y_vector(2,:)-r(2,:)).^2)+sqrt((y_vector(3,:)-r(3,:)).^2)+sqrt((y_vector(4,:)-r(4,:)).^2))
% ����ʱ��
elapsed_time=cputime-t0