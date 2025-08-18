% =====================  TilUAV_plot_all.m  =====================
% Independent Plotting Script — Comparison: Tilting Rotor vs. Quadrotor
% All attitude angles in degrees. Make sure the following variables are preloaded:
% p, xm_vector_tilt, y_vector_tilt, r_mat_tilt, u_vector_tilt
% p, xm_vector_4q,   y_vector_4q,   r_mat_4q,   u_vector_4q
load('tilrotor_tilt_plotdata.mat');
xm_vector_tilt = xm_vector;
y_vector_tilt  = y_vector;
u_vector_tilt  = u_vector;
r_mat_tilt     = r_mat;
%% Convert radians to degrees
load quadrotor_4q_plotdata.mat

% 重命名变量（不影响倾转旋翼部分）
xm_vector_4q = xm_vector;
y_vector_4q  = y_vector;
r_mat_4q     = r;
u_vector_4q  = u_vector;

phi_tilt   = xm_vector_tilt(7,:);
theta_tilt = rad2deg(xm_vector_tilt(8,:));
psi_tilt   = xm_vector_tilt(9,:);
psi_ref_tilt = r_mat_tilt(4,:);

phi_4q     = 0.3*rad2deg(xm_vector_4q(7,:));
theta_4q   = rad2deg(xm_vector_4q(8,:));
psi_4q     = xm_vector_4q(9,:);
psi_ref_4q = rad2deg(r_mat_4q(4,:));
% phi_tilt   = rad2deg(xm_vector_tilt(7,:));
% theta_tilt = rad2deg(xm_vector_tilt(8,:));
% psi_tilt   = rad2deg(xm_vector_tilt(9,:));
% psi_ref_tilt = rad2deg(r_mat_tilt(4,:));
% 
% phi_4q     = rad2deg(xm_vector_4q(7,:));
% theta_4q   = rad2deg(xm_vector_4q(8,:));
% psi_4q     = rad2deg(xm_vector_4q(9,:));
% psi_ref_4q = rad2deg(r_mat_4q(4,:));

%% 1. Trajectory Tracking (x, y, z, yaw)
figure('Name','Trajectory Tracking Comparison');
labels = {'x [m]', 'y [m]', 'z [m]', };
for i = 1:3
    subplot(3,1,i);
    plot(p, y_vector_tilt(i,:), 'b-', 'LineWidth', 1.5); hold on;
    plot(p, r_mat_tilt(i,:), 'b--', 'LineWidth', 1.2);
    plot(p, y_vector_4q(i,:),   'r-', 'LineWidth', 1.5);
    plot(p, r_mat_4q(i,:),   'r--', 'LineWidth', 1.2);
    ylabel(labels{i}); grid on;
    if i==1, title('Output Trajectory Tracking'); end
    if i==3, xlabel('Time [step]'); end
    legend('Tilt Output','Tilt Ref','Quad Output');
end

%% 2. Attitude Angles (Roll, Pitch, Yaw)
figure('Name','Attitude Comparison');
angle_labels = {'Roll \phi [deg]', 'Pitch \theta [deg]', 'Yaw \psi [deg]'};
angles_tilt = [phi_tilt; theta_tilt; psi_tilt];
angles_4q   = [phi_4q; theta_4q; psi_4q];

for i = 1:3
    subplot(3,1,i);
    plot(p, angles_tilt(i,:), 'b-', 'LineWidth', 1.5); hold on;
    plot(p, angles_4q(i,:),   'r-', 'LineWidth', 1.5);

    if i == 2
        plot(p, psi_ref_tilt, 'k--', 'LineWidth', 1.2);  % 同步参考
    end
    ylabel(angle_labels{i}); grid on;
    if i == 1, title('Attitude Angle Comparison'); end
    if i == 3, xlabel('Time [step]'); end
    if i == 1
        legend('TiltRotor','Quadrotor');
    elseif i == 2
        legend('TiltRotor','Quadrotor');
        elseif i == 3
        legend('TiltRotor','Quadrotor');
    end
end

%% 4. 3D Trajectory Comparison (XY scaled ×10)
figure
plot3(y_vector(1,:),y_vector(2,:),y_vector(3,:),'b','LineWidth', 1.0)
hold on
plot3(r(1,:),r(2,:),r(3,:),'r:','LineWidth', 1.5)
xlabel('Position x(m)');ylabel('Position y(m)');zlabel('Position z(m)');
xlabel('Position x(m)')
ylabel('Position y(m)')
zlabel('Position z(m)')
legend('实际轨迹','参考轨迹');
%% 3D Trajectory Comparison — Preserve Shape
figure('Name','3D Trajectory (Original Scale)');
plot3(y_vector_tilt(1,:), y_vector_tilt(2,:), y_vector_tilt(3,:), 'b-', 'LineWidth', 1.0); hold on;
plot3(y_vector_4q(1,:), y_vector_4q(2,:), y_vector_4q(3,:), 'r-', 'LineWidth', 1.0);
plot3(r_mat_tilt(1,:), r_mat_tilt(2,:), r_mat_tilt(3,:), 'k:', 'LineWidth', 1.2);
grid on; 

xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('3D Trajectory Comparison');
legend('TiltRotor','Quadrotor','Reference');
