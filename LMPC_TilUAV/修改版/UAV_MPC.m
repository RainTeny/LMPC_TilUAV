
%% ���������켣���������MPC��ģ�����
clear all
close all
t0=cputime;
%% ��һ���֣�Quadcopter state-space model
%% ϵͳ����
m=2; % mass kg
g=9.81; 
Ix=7.5E-3; % ת������ x �� kg*m^2
Iy=7.5E-3; % ת������ y �� kg*m^2
Iz=1.3E-2; % ת������ z �� kg*m^2
dx=0; dy=0; dz=0; % x y z �����Ŷ�Ϊ��
kx=0.25; ky=0.25; kz=0.25; % ��������

%% ״̬�ռ�ģ��
%δ�򻯵�ģ��
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12
syms u1 u2 u3 u4
% x=[x y z vx vy vz roll pitch yaw wr wp wy]
xdot1=x4;
xdot2=x5;
xdot3=x6;
xdot4=(u1/m)*(cos(x9)*sin(x8)*cos(x7)+sin(x9)*sin(x7))-kx*x4/m+dx/m;
xdot5=(u1/m)*(sin(x9)*sin(x8)*cos(x7)-cos(x9)*sin(x7))-ky*x5/m+dy/m;
xdot6=(u1/m)*cos(x8)*cos(x7)-g-kz*x6/m+dz/m;
xdot7=x10+x11*sin(x7)*tan(x8)+x12*cos(x7)*tan(x8);
xdot8=x11*cos(x7)-x12*sin(x7);
xdot9=sin(x7)*x11/cos(x8)+cos(x7)*x12/cos(x8);
xdot10=(u2/Ix)-((Iy-Iz)/Ix)*x11*x12;
xdot11=(u3/Iy)-((Iz-Ix)/Iy)*x10*x12;
xdot12=(u4/Iz)-((Ix-Iy)/Iz)*x10*x11;
xdot=[xdot1,xdot2,xdot3,xdot4,xdot5,xdot6,xdot7,xdot8,xdot9,xdot10,xdot11,xdot12]';
x=[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12]';
u=[u1,u2,u3,u4]';
y=[x1,x2,x3,x7,x8,x9]';

%% �򻯺��ģ��
% ���ƿ��� sin x -> x; cos x -> 1; tan x -> x; no distrurbances
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12
syms u1 u2 u3 u4
xdot1=x4;
xdot2=x5;
xdot3=x6;
xdot4=(u1/m)*(x8+x9*x7)-kx*x4/m;
xdot5=(u1/m)*(x9*x8-x7)-ky*x5/m;
xdot6=(u1/m)-g-kz*x6/m;
xdot7=x10+x11*x7*x8+x12*x8;
xdot8=x11-x12*x7;
xdot9=x7*x11+x12;
xdot10=(u2/Ix)-((Iy-Iz)/Ix)*x11*x12;
xdot11=(u3/Iy)-((Iz-Ix)/Iy)*x10*x12;
xdot12=(u4/Iz)-((Ix-Iy)/Iz)*x10*x11;
xdot=[xdot1 xdot2 xdot3 xdot4 xdot5 xdot6 xdot7 xdot8 xdot9 xdot10 xdot11 xdot12].';
x=[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12].';
u=[u1,u2,u3,u4].';
y=[x1,x2,x3,x9].';
[x_size,aux]=size(x);
[y_size,aux]=size(y);
[u_size,aux]=size(u);

%% �ſɱ����Ի�
% �ȶ���
ue=[m*g,0,0,0].';
xe=[x1,x2,x3,0,0,0,0,0,0,0,0,0].';%״̬����ʼֵ
JA=jacobian(xdot,x.');
JB=jacobian(xdot,u.');
JC=jacobian(y, x.');
A=subs(JA,[x;u],[xe;ue]); A=eval(A);%�����ȶ���
B=subs(JB,[x;u],[xe;ue]); B=eval(B);
C=subs(JC,[x;u],[xe;ue]); C=eval(C);
%%��ɢ״̬�ռ�ģ��
fs=50; % 100 Hz
Ts=1/fs; 
sysc=ss(A,B,C,0);  %%����ʵ��������״̬�ռ�ģ��
sysd=c2d(sysc,Ts); %%��s��ı��ʽת����z��ı��ʽ
Am=sysd.A;
Bm=sysd.B;
Cm=sysd.C;
Dm=zeros(x_size,1);

%% MPC������
%% ��ʼ��
N_sim=10*fs; %% samples = seconds*frequency
Nc=10;  % ����ʱ��
Np=50; % Ԥ��ʱ��
R = 0.01;   %������Ȩ��

% Init control, reference and output signal
u=zeros(u_size,1);  
y=zeros(y_size,1); 

rx=0*ones(1,N_sim);
ry=0*ones(1,N_sim);
rz=0*ones(1,N_sim);
ryaw=0*ones(1,N_sim);

xm_vector=[];
u_vector=[];
deltau_vector=[];
y_vector=[];

% ��ʼϵͳ״̬ ȫ��Ϊ0
xm=zeros(x_size,1); % ״̬����
Xf=zeros(x_size+y_size,1); % ��ǿ����״̬ [deltax y]'

%% �õ���������ģ�ͺ��������Ʋ���
[GG, GF, GR, A, B, C , GD, F, G]=mpcgain_mimo1(Am,Bm,Cm,Nc,Np,Dm); 

%% �����켣���Ƶĺ㶨����
f1=GG+R*eye(Nc*u_size,Nc*u_size); % E for the cost function J=xEx'+x'F

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
M_deltau=[aux(1,:);-aux(1,:);
    aux(2,:);-aux(2,:);
    aux(3,:);-aux(3,:);
    aux(4,:);-aux(4,:)];

M_u=[aux(1,:);-aux(1,:);
    aux(2,:);-aux(2,:);
    aux(3,:);-aux(3,:);
    aux(4,:);-aux(4,:)];

M=[M_output; M_deltau; M_u];

%% ѭ��
for k=1:N_sim 
    %% �ο��켣
    rx(k)=sin(0.03*k);    
    ry(k)=cos(0.03*k);       
    rz(k)=0.1*k;    
    ryaw(k)=30;
    
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
    deltau=DeltaU(1:4);  % Ӧ�ù����Ż�����
    u=u+deltau; % (u(k)=u(k-1)+deltau(k))
    u=u-[m*g;0;0;0]; % �ȶ���
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
  
end
% ������ͨ�����������������ڶԱȻ�ͼ��
save('quadrotor_4q_data.mat', 'xm_vector', 'y_vector', 'r', 'u_vector', 'N_sim');

% ����Ϊͳһ����
save('quadrotor_4q_plotdata.mat', ...
    'xm_vector', 'y_vector', 'u_vector', ...
    'r', 'N_sim');



UAV_plot;

% ���ڼ����ۼ�ʱ�����
errorsum=sum(sqrt((y_vector(1,:)-r(1,:)).^2)+sqrt((y_vector(2,:)-r(2,:)).^2)+sqrt((y_vector(3,:)-r(3,:)).^2)+sqrt((y_vector(4,:)-r(4,:)).^2))
% ����ʱ��
elapsed_time=cputime-t0