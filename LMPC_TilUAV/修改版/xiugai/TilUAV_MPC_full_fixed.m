
% =====================  TilUAV_MPC_full_fixed.m  =====================
% 修复版本：线性增量式MPC控制器 + 俯仰角跟踪，含姿态图与注释
clear; close all; clc;
t0 = cputime;

%% 1  系统参数 -----------------------------------------------------------
m  = 1;              lm = 0.18;           d  = 1.574e-7;  g  = 9.81;
Ix = 7.5e-3;  Iy = 7.5e-3;  Iz = 1.3e-2;   kx = .25;  ky = .25;  kz = .25;

%% 2  简化模型 ----------------------------------------------------------
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 u1 u2 u3 u4 u5 u6
xdot = [x4;
         x5;
         x6;
         (u1+u2)*cos(x8)*cos(x9) + (u3+u4+u5+u6)*(cos(x7)*sin(x8)*cos(x9)+sin(x7)*sin(x9)) - kx*x4;
         (u1+u2)*cos(x8)*sin(x9) + (u3+u4+u5+u6)*(cos(x7)*sin(x8)*cos(x9)-sin(x7)*cos(x9)) - ky*x5;
        -(u1+u2)*sin(x8) + (u3+u4+u5+u6)*cos(x7)*cos(x8) - g - kz*x6;
         x10 + x11*sin(x7)*x8 + x12*cos(x7)*x8;
         x11*cos(x7) - x12*sin(x7);
         sin(x7)*x11 + cos(x7)*x12;
         (lm*(u3-u4)/Ix) - ((Iy-Iz)/Ix)*x11*x12;
         (lm*(u5-u6)/Iy) - ((Iz-Ix)/Iy)*x10*x12;
         (lm*(u1+u2)/Iz) + (d*(u3+u4-u5-u6)/Iz) - ((Ix-Iy)/Iz)*x10*x11];

x = [x1;x2;x3;x4;x5;x6;x7;x8;x9;x10;x11;x12];
u = [u1;u2;u3;u4;u5;u6];
y_sym = [x1;x2;x3;x8];

x_size = numel(x); y_size = numel(y_sym); u_size = numel(u);

%% 3  线性化与离散 -------------------------------------------------------
ue = [0;0;0.25*g;0.25*g;0.25*g;0.25*g];
A = double(subs(jacobian(xdot,x), [x;u], [zeros(x_size,1);ue]));
B = double(subs(jacobian(xdot,u), [x;u], [zeros(x_size,1);ue]));
C = double(subs(jacobian(y_sym,x), [x;u], [zeros(x_size,1);ue]));
fs = 50; Ts = 1/fs; sysd = c2d(ss(A,B,C,0),Ts);
Am = sysd.A; Bm = sysd.B; Cm = sysd.C;
Dm = zeros(size(Cm,1),1);

%% 4  MPC 设置 -------------------------------------------------------------
Nc = 15; Np = 25;
Wy = diag([100 100 10 100]);
Wu = diag([0.2 0.2 0.1 0.1 0.1 0.1]);
Qy = kron(eye(Np), Wy);
R = kron(eye(Nc), Wu);

[Phi_Phi, Phi_F, Phi_R, ~, ~, ~, ~, ~, ~] = mpcgain_mimo(Am,Bm,Cm,Nc,Np,Dm,Qy);

%% 5  仿真准备 -----------------------------------------------------------
N_sim = 10*fs; p = 0:N_sim-1;
rx = zeros(1,N_sim); ry = rx; rz = rx; theta_ref = rx;

xm = zeros(x_size,1); Xf = zeros(x_size+y_size,1);
xm_vector = zeros(x_size,N_sim);
y_vector = zeros(y_size,N_sim);
u_vector = zeros(u_size,N_sim);
deltau_vector = zeros(u_size,N_sim);

%% 6  主循环 -------------------------------------------------------------
for k = 1:N_sim
    rx(k) = sin(0.03*k); ry(k) = cos(0.03*k); rz(k) = 0.1*k*Ts*fs;
    theta_ref(k) = 0;  % 恒定期望俯仰角

    r_k = [rx(k); ry(k); rz(k); theta_ref(k)];

    f1 = Phi_Phi + R;
    f2 = Phi_R*r_k - Phi_F*Xf;
    DeltaU = f1 \ f2;
    deltau = DeltaU(1:u_size);

    if k == 1, u = zeros(u_size,1); end
    u = u + deltau; u = u - ue;
    xm_old = xm;
    xm = Am*xm + Bm*u;
    y = Cm*xm;
    Xf = [xm - xm_old; y];

    xm_vector(:,k) = xm;
    y_vector(:,k) = y;
    u_vector(:,k) = u;
    deltau_vector(:,k) = deltau;
end

r_mat = [rx; ry; rz; theta_ref];
err = sum(sqrt((y_vector(1,:)-rx).^2) + sqrt((y_vector(2,:)-ry).^2) + ...
          sqrt((y_vector(3,:)-rz).^2) + sqrt((y_vector(4,:)-theta_ref).^2));


%% 7  绘图输出 -------------------------------------------------------------
figure;
for i = 1:4
    subplot(4,1,i);
    plot(p, y_vector(i,:), 'b', p, r_mat(i,:), 'r--', 'LineWidth', 1);
    grid on;
end
subplot(4,1,1); legend('x','x_d'); subplot(4,1,2); legend('y','y_d');
subplot(4,1,3); legend('z','z_d'); subplot(4,1,4); legend('	heta','	heta_d');

figure; hold on;
plot(p, xm_vector(7,:), 'r-', 'LineWidth', 1.5);  % roll
plot(p, xm_vector(8,:), 'g-', 'LineWidth', 1.5);  % pitch
plot(p, xm_vector(9,:), 'b-', 'LineWidth', 1.5);  % yaw
legend('\phi (roll)', '	heta (pitch)', '\psi (yaw)');
xlabel('时间步'); ylabel('角度(rad)'); title('姿态角'); grid on;
