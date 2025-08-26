function [X, U, X1, X2, U1, U2] = get_rnd_trajectories(X0, n_control, t_traj, show_plot, flag)
% 生成随机短轨迹（两段常值控制），围绕悬停点采样，
% 并覆盖上升/下降/前飞/侧飞等不同速度段
% 过滤规则：
%   - 竖直合力 Fz >= hover_margin * m*g（两段都要满足）
%   - 姿态正常：|phi|、|theta| < 20 deg（整个轨迹；硬限）
% 输出矩阵按样本严格对齐

%% 参数
params = get_params();
m  = params.mass;
g  = params.g;
hover_thrust = m*g;                         % 悬停总推力

T   = numel(t_traj);
Tm1 = T-1;
mid = max(2, floor(T/2));                   % 两段分界（至少保证>1）

% —— 轻度收敛随机度 + 小姿态初值（更便于可辨识）
switch flag
    case 'train'
        sigma_total  = 0.40;                 % 原 0.20 -> 0.12
        sigma_split  = 0.80;                 % 原 0.50 -> 0.40
        alpha_sigma1 = deg2rad(30);          % 原 20° -> 14°
        alpha_sigma2 = deg2rad(30);
        hover_margin = 0.88;
        v0_horz_max  = 0.3;
        v0_vert_max  = 0.5;
        att0_max_deg = 3;                   % 原 10 -> 12（更接近真实微姿态）
    case 'val'
        sigma_total  = 0.40;
        sigma_split  = 0.80;
        alpha_sigma1 = deg2rad(30);
        alpha_sigma2 = deg2rad(30);
        hover_margin = 0.88;
        v0_horz_max  = 0.3;
        v0_vert_max  = 0.5;
        att0_max_deg = 3;
    otherwise  % 'mpc'：更温和的扰动
        sigma_total  = 0.08;
        sigma_split  = 0.30;
        alpha_sigma1 = deg2rad(12);
        alpha_sigma2 = deg2rad(12);
        hover_margin = 0.98;
        v0_horz_max  = 0.50;
        v0_vert_max  = 0.30;
        att0_max_deg = 6;
end

% 硬限制（为防止过度注入）
alpha_lim = deg2rad(20);     % 倾转角硬限幅（±20°）
ang_lim   = deg2rad(20);     % |phi|、|theta| < 20°

% “定向激励”注入概率与强度
p_yaw_pulse   = 0.6;               % 60% 轨迹加一个 yaw 轻量脉冲（同相）
p_roll_pulse  = 0.5;               % 50% 轨迹加一个 roll 轻量脉冲（反相）
p_pitch_pulse = 0.5;               % 50% 轨迹加一个 pitch 轻量脉冲（F3/F4 差分）
yaw_amp   = deg2rad(6);            % ≈6°
roll_amp  = deg2rad(5);            % ≈5°
dF_pitch  = 0.06 * (m*g);          % ≈6%·mg 小差分

%% 容器
X  = []; X1 = []; X2 = [];
U  = []; U1 = []; U2 = [];
keep  = 0;  bad_theta = 0; bad_Fz = 0;

% 预先拆分时间轴（两段）
t1 = t_traj(1:mid);
t2 = t_traj(mid:end);

%% 循环采样
for i = 1:n_control

    % -------- 随机初始状态（速度+小姿态） --------
    x0 = X0;
    v0_xy = v0_horz_max * (2*rand(2,1)-1);      % [-v0_horz_max, v0_horz_max]
    v0_z  = v0_vert_max * (2*rand-1);
    x0(4:6) = [v0_xy; v0_z];

    att0 = deg2rad(att0_max_deg) * (2*rand(3,1)-1); % [phi0, theta0, psi0]
    % 限制 |phi|、|theta|
    att0(1) = clip(att0(1), -ang_lim, ang_lim); % phi
    att0(2) = clip(att0(2), -ang_lim, ang_lim); % theta
    x0(7:9) = att0;

    % -------- 设计两段常值控制，制造不同速度段 --------
    % 0: 悬停微扰 → 前飞
    % 1: 前飞 → 减速
    % 2: 上升 → 前飞
    % 3: 侧飞 → 回悬停
    % 4: 下降 → 侧飞
    mode = randi([0,4],1);

    % —— 段1：总推力
    thrust1 = hover_thrust * (1 + sigma_total*randn);
    base1   = thrust1/4;
    eps1    = sigma_split * base1 * randn(4,1);
    eps1    = eps1 - mean(eps1);
    F1      = max(base1 + eps1, 0);

    % —— 段2：总推力
    thrust2 = hover_thrust * (1 + sigma_total*randn);
    base2   = thrust2/4;
    eps2    = sigma_split * base2 * randn(4,1);
    eps2    = eps2 - mean(eps2);
    F2      = max(base2 + eps2, 0);

    % 倾转角（两个可倾转桨：1、3）
    switch mode
        case 0  % 悬停微扰 → 前飞（a1=a2 同向前）
            a1_1 = clip(alpha_sigma1*randn, -alpha_lim, alpha_lim)*0.3;
            a2_1 = clip(alpha_sigma1*randn, -alpha_lim, alpha_lim)*0.3;
            a1_2 = clip(alpha_sigma2*randn + deg2rad(8), -alpha_lim, alpha_lim);
            a2_2 = clip(alpha_sigma2*randn + deg2rad(8), -alpha_lim, alpha_lim);

        case 1  % 前飞 → 减速（先前倾，再回零附近）
            a1_1 = clip(alpha_sigma1*randn + deg2rad(8), -alpha_lim, alpha_lim);
            a2_1 = clip(alpha_sigma1*randn + deg2rad(8), -alpha_lim, alpha_lim);
            a1_2 = clip(alpha_sigma2*randn, -alpha_lim, alpha_lim)*0.2;
            a2_2 = clip(alpha_sigma2*randn, -alpha_lim, alpha_lim)*0.2;

        case 2  % 上升 → 前飞（先 thrust↑，再前倾）
            thrust1 = max(thrust1,  hover_thrust*(1 + 0.10 + 0.05*rand));  % 上升
            base1   = thrust1/4; F1 = max(base1 + eps1, 0);
            a1_1 = clip(alpha_sigma1*randn, -alpha_lim, alpha_lim)*0.2;
            a2_1 = clip(alpha_sigma1*randn, -alpha_lim, alpha_lim)*0.2;
            a1_2 = clip(alpha_sigma2*randn + deg2rad(8), -alpha_lim, alpha_lim);
            a2_2 = clip(alpha_sigma2*randn + deg2rad(8), -alpha_lim, alpha_lim);

        case 3  % 侧飞 → 回悬停（a1=-a2 产生侧向分量）
            sgn   = sign(randn); if sgn==0, sgn=1; end
            a1_1 = clip(alpha_sigma1*randn + sgn*deg2rad(8), -alpha_lim, alpha_lim);
            a2_1 = clip(alpha_sigma1*randn - sgn*deg2rad(8), -alpha_lim, alpha_lim);
            a1_2 = clip(alpha_sigma2*randn, -alpha_lim, alpha_lim)*0.3;
            a2_2 = clip(alpha_sigma2*randn, -alpha_lim, alpha_lim)*0.3;

        otherwise % 4: 下降 → 侧飞
            thrust1 = min(thrust1,  hover_thrust*(1 - 0.08*rand));         % 下降
            base1   = thrust1/4; F1 = max(base1 + eps1, 0);
            sgn   = sign(randn); if sgn==0, sgn=1; end
            a1_1 = clip(alpha_sigma1*randn, -alpha_lim, alpha_lim)*0.2;
            a2_1 = clip(alpha_sigma1*randn, -alpha_lim, alpha_lim)*0.2;
            a1_2 = clip(alpha_sigma2*randn + sgn*deg2rad(8), -alpha_lim, alpha_lim);
            a2_2 = clip(alpha_sigma2*randn - sgn*deg2rad(8), -alpha_lim, alpha_lim);
    end

    % ====== 在两段上注入“定向转动激励”（常值脉冲，仍保持两段常值）======
    % yaw 轻量脉冲（同相，激励 Mz）
    if rand < p_yaw_pulse
        if rand < 0.5
            a1_1 = clip(a1_1 + yaw_amp*(2*(rand>0.5)-1), -alpha_lim, alpha_lim);
            a2_1 = clip(a2_1 + yaw_amp*(2*(rand>0.5)-1), -alpha_lim, alpha_lim);
        else
            a1_2 = clip(a1_2 + yaw_amp*(2*(rand>0.5)-1), -alpha_lim, alpha_lim);
            a2_2 = clip(a2_2 + yaw_amp*(2*(rand>0.5)-1), -alpha_lim, alpha_lim);
        end
    end

    % roll 轻量脉冲（反相，激励 Mx）
    if rand < p_roll_pulse
        if rand < 0.5
            a1_1 = clip(a1_1 +  roll_amp, -alpha_lim, alpha_lim);
            a2_1 = clip(a2_1 -  roll_amp, -alpha_lim, alpha_lim);
        else
            a1_2 = clip(a1_2 +  roll_amp, -alpha_lim, alpha_lim);
            a2_2 = clip(a2_2 -  roll_amp, -alpha_lim, alpha_lim);
        end
    end

    % pitch 轻量脉冲（F3/F4 差分，激励 My）
    if rand < p_pitch_pulse
        if rand < 0.5
            F1 = F1; % no change
            F2 = F2; % no change
            F3 = max(F1(3) + dF_pitch, 0);
            F4 = max(F1(4) - dF_pitch, 0);
            F1 = [F1(1); F1(2); F3; F4];
        else
            F3 = max(F2(3) + dF_pitch, 0);
            F4 = max(F2(4) - dF_pitch, 0);
            F2 = [F2(1); F2(2); F3; F4];
        end
    end

    % 两段竖直分量校验（注入后再校验）
    Fz1 = cos(a1_1)*F1(1) + F1(2) + cos(a2_1)*F1(3) + F1(4);
    Fz2 = cos(a1_2)*F2(1) + F2(2) + cos(a2_2)*F2(3) + F2(4);
    if (Fz1 < hover_margin*hover_thrust) || (Fz2 < hover_margin*hover_thrust)
        bad_Fz = bad_Fz + 1; 
        continue;
    end

    Ui1 = [F1; a1_1; a2_1];
    Ui2 = [F2; a1_2; a2_2];

    % -------- 两段仿真并拼接（仍是两次 ode45） --------
    [~, x1] = ode45(@(t,X) dynamics_SRB(t, X, Ui1, params), t1, x0);
    x1e = x1(end,:).';
    [~, x2] = ode45(@(t,X) dynamics_SRB(t, X, Ui2, params), t2, x1e);
    x = [x1; x2(2:end,:)];

    % -------- 姿态过滤：|phi|、|theta| < 20° --------
    phi_seq   = x(:,7);
    theta_seq = x(:,8);
    if any(abs(phi_seq) > ang_lim) || any(abs(theta_seq) > ang_lim)
        bad_theta = bad_theta + 1;
        continue;
    end

    % -------- 通过筛选，入库 --------
    keep = keep + 1;

    % 状态堆叠（列拼）
    X  = [X,  x'];                 % 12×T
    X1 = [X1, x(1:end-1,:)'];      % 12×(T-1)
    X2 = [X2, x(2:end,:)'];        % 12×(T-1)

    % 控制按时间对齐（两段常值）
    u_seq_full  = [ repmat(Ui1,1,mid), repmat(Ui2,1,T-mid) ];
    u_seq_short = [ repmat(Ui1,1,mid-1), repmat(Ui2,1,(T-mid)) ];

    U  = [U,  u_seq_full      ];   % 6×T
    U1 = [U1, u_seq_short     ];   % 6×(T-1)
    U2 = [U2, u_seq_short     ];

    % 可视化
    if show_plot
        if keep==1
            figure(99); clf; hold on; grid on; axis equal; view(3);
            xlabel('x'); ylabel('y'); zlabel('z');
            title('EDMD Training Trajectories (tilt-rotor, 2-phase inputs + light pulses)');
        end
        plot3(x(:,1), x(:,2), x(:,3), 'b-');
        drawnow limitrate
    end
end

%% 日志
fprintf('\n筛选结果：\n');
fprintf('  通过轨迹：%d / %d\n', keep, n_control);
fprintf('  剔除(姿态超限)：%d  |  剔除(Fz不足)：%d  | 阈值 Fz >= %.2f mg\n', ...
        bad_theta, bad_Fz, hover_margin);

if keep==0
    warning('没有任何轨迹通过筛选，返回空矩阵。');
end
end

% --- 小工具函数 ---
function y = clip(x, lo, hi)
    y = min(max(x, lo), hi);
end

