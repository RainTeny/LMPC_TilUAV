function dXdt = dynamics_SRB(t, X, U, params)
% Modified: EDMD-compatible Euler-based Tilt Rotor dynamics
%% Parameters
m  = params.mass;
g  = 9.81;
lm = params.lm;
Ix = params.J(1,1);
Iy = params.J(2,2);
Iz = params.J(3,3);
kx = params.kx;
ky = params.ky;
kz = params.kz;
d  = params.d;

%% State unpack
x1 = X(1); x2 = X(2); x3 = X(3);
dx = X(4); dy = X(5); dz = X(6);
phi = X(7); theta = X(8); psi = X(9);
p = X(10); q = X(11); r = X(12);


%% Input unpack
F1 = U(1); F2 = U(2); F3 = U(3); F4 = U(4);
alpha1 = U(5); alpha2 = U(6);

%% Tilt rotor force decomposition
F_tilt_xy = sin(alpha1)*F1 + sin(alpha2)*F3;
F_total_z = cos(alpha1)*F1 + cos(alpha2)*F3 + F2 + F4;

%% Translational dynamics
xdot1 = dx;
xdot2 = dy;
xdot3 = dz;
% 使用机体系中的推力方向（Z轴向上）作用力
F_body = [0; 0; F_total_z];

% 构建旋转矩阵（ZYX顺序，即 ψ→θ→φ）
R = eul2rotm([psi, theta, phi], 'ZYX');

% 变换为惯性坐标系
F_inertial = R * F_body;

% 平动加速度：F=ma 减阻力
xdot4 = F_inertial(1)/m - kx*dx/m;
xdot5 = F_inertial(2)/m - ky*dy/m;
xdot6 = F_inertial(3)/m - g - kz*dz/m;


%% Rotational kinematics (Euler)
% xdot7 = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
% xdot8 = q*cos(phi) - r*sin(phi);
% xdot9 = q*sin(phi)/(cos(theta)+1e-6) + r*cos(phi)/(cos(theta)+1e-6);
xdot7 = p ;
xdot8 = q;
xdot9 = r;
%% Angular dynamics
xdot10 = (lm*(cos(alpha1)*F1 - cos(alpha2)*F3))/Ix - ((Iy - Iz)/Ix)*q*r;
xdot11 = (lm*(F3 - F4))/Iy - ((Iz - Ix)/Iy)*p*r;
xdot12 = (lm*(sin(alpha1)*F1 + sin(alpha2)*F3))/Iz + ...
         (d*(cos(alpha1)*F1 + cos(alpha2)*F3 - F3 - F4))/Iz - ...
         ((Ix - Iy)/Iz)*p*q;
%%检查错误    
% if F_total_z < m*g
%     disp(['Fz不够:', num2str(F_total_z), ' < ', num2str(m*g)]);
% end


%% Return
dXdt = [xdot1; xdot2; xdot3;
        xdot4; xdot5; xdot6;
        xdot7; xdot8; xdot9;
        xdot10; xdot11; xdot12];
end








