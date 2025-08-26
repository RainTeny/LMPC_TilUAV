function STATS = eval_EDMD(X0, dt, t_span, EDMD, show_plot, seed)
% 评估训练好的 EDMD 模型（与训练一致）
% 保留：RMSE + 论文口径 nRMSE%
% 新增：多步误差曲线 + R²
% 去掉：Ljung–Box 自相关检验、谱半径/特征值可视化

if nargin < 5 || isempty(show_plot), show_plot = false; end
if nargin < 6 || isempty(seed)
    rng('shuffle');
else
    rng(seed);
end

% ==== 模型参数 ====
A = EDMD.A;  B = EDMD.B;  C = EDMD.C;
E = [];  if isfield(EDMD,'E')  && ~isempty(EDMD.E),  E  = EDMD.E;  end
Bj = {}; if isfield(EDMD,'Bj') && ~isempty(EDMD.Bj), Bj = EDMD.Bj; end
b = zeros(size(A,1),1); if isfield(EDMD,'b') && ~isempty(EDMD.b), b = EDMD.b; end

use_cross       = ~isempty(E) || ~isempty(Bj);
cross_with_phys = isfield(EDMD,'opts') && isfield(EDMD.opts,'cross_with_phys') && EDMD.opts.cross_with_phys;
strip_const     = isfield(EDMD,'opts') && isfield(EDMD.opts,'strip_const')     && EDMD.opts.strip_const;

% ==== 时间轴 ====
t_grid = 0:dt:t_span;
T = numel(t_grid);
Nstep = T - 1;

% ==== 评估采样 ====
pool_try  = 1000;
K_eval    = 50;                      % 评估轨迹条数
theta_thr = deg2rad(87);
[ X_raw, U_raw, ~, ~, ~, ~ ] = get_rnd_trajectories(X0, pool_try, t_grid, false, 'val');

X_list = {}; U_list = {};
for k = 1:pool_try
    if numel(X_list) >= K_eval, break; end
    cols = (k-1)*T + 1 : k*T;
    if all(abs(X_raw(8, cols)) < theta_thr)
        X_list{end+1} = X_raw(:, cols);
        U_list{end+1} = U_raw(:, cols(1:end-1));
    end
end
assert(~isempty(X_list), '评估采样无合格轨迹，请放宽阈值或增大采样池。');
fprintf('✔ 评估：采到 %d 条合格轨迹\n', numel(X_list));

% ==== 画图准备 ====
if show_plot
    figure(11); clf; hold on; grid on; axis tight; view(3);
    xlabel('x'); ylabel('y'); zlabel('z');
    title('True (blue) vs Predicted (red)');
end
Mplot = 50;

% ==== 构造 z 的工具 ====
nZ   = size(C,2);
nPhi = nZ - 12;
make_z = @(x) local_make_z_guard(x, nPhi, strip_const);

% ==== 统计容器 ====
rmse_pos  = zeros(1, numel(X_list));
rmse_full = zeros(1, numel(X_list));

num_all = struct('p',0,'v',0,'t',0,'w',0);
den_all = struct('p',0,'v',0,'t',0,'w',0);
traj_pct = struct('p',[],'v',[],'t',[],'w',[]);


% R² 收集
pred_store = cell(1, numel(X_list));

for n = 1:numel(X_list)
    X_true = X_list{n};   % 12×T
    U_seq  = U_list{n};   %  6×(T-1)

    % 初值
    x = X_true(:,1);
    z = make_z(x);

    if n == 1
        fprintf('[诊断] nZ=%d, nPhi=%d, size(C)=%dx%d\n', nZ, nPhi, size(C,1), size(C,2));
        if use_cross && ~isempty(E)
            if cross_with_phys, zu_probe = build_cross(z, U_seq(:,1), x);
            else,               zu_probe = build_cross(z, U_seq(:,1), []);
            end
            assert(size(E,2)==numel(zu_probe), 'E*zu 维度不匹配：size(E,2)=%d, numel(zu)=%d。', size(E,2), numel(zu_probe));
        end
    end

    X_pred = zeros(12, T);
    X_pred(:,1) = x;

    % === 时间推进 ===
    for i = 1:Nstep
        u = U_seq(:,i);
        x_cur = C*z;

        if use_cross
            if ~isempty(Bj)     % 双线性分块
                uz_term = zeros(nZ,1);
                for j = 1:numel(Bj)
                    uz_term = uz_term + u(j) * (Bj{j} * z);
                end
            else                % 退回 E*(u⊗z)
                if cross_with_phys, zu = build_cross(z, u, x_cur);
                else,               zu = build_cross(z, u, []);
                end
                uz_term = E*zu;
            end
            z_next = A*z + B*u + uz_term + b;
        else
            z_next = A*z + B*u + b;
        end

        if ~all(isfinite(z_next)) || norm(z_next) > 1e8
            X_pred(:,i+1:end) = repmat(C*z, 1, T-i);  % 兜底
            break;
        end

        x_next = C*z_next;
        X_pred(:,i+1) = x_next;
        z = make_z(x_next);


    end

    % 误差
    err = X_true - X_pred;
    err_pos  = err(1:3,:);
    rmse_pos(n)  = sqrt(mean(err_pos(:).^2));
    rmse_full(n) = sqrt(mean(err(:).^2));

    pred_store{n} = X_pred;

    % 论文口径 nRMSE% 分组
    idx.p = 1:3;  idx.v = 4:6;  idx.t = 7:9;  idx.w = 10:12;
    num_traj = structfun(@(I) sum(vecnorm(err(I,:),2,1).^2),    idx, 'UniformOutput', false);
    den_traj = structfun(@(I) sum(vecnorm(X_true(I,:),2,1).^2), idx, 'UniformOutput', false);

    traj_pct.p(end+1) = 100 * sqrt( safe_div(num_traj.p, den_traj.p) ); %#ok<AGROW>
    traj_pct.v(end+1) = 100 * sqrt( safe_div(num_traj.v, den_traj.v) ); %#ok<AGROW>
    traj_pct.t(end+1) = 100 * sqrt( safe_div(num_traj.t, den_traj.t) ); %#ok<AGROW>
    traj_pct.w(end+1) = 100 * sqrt( safe_div(num_traj.w, den_traj.w) ); %#ok<AGROW>

    num_all.p = num_all.p + num_traj.p;  den_all.p = den_all.p + den_traj.p;
    num_all.v = num_all.v + num_traj.v;  den_all.v = den_all.v + den_traj.v;
    num_all.t = num_all.t + num_traj.t;  den_all.t = den_all.t + den_traj.t;
    num_all.w = num_all.w + num_traj.w;  den_all.w = den_all.w + den_traj.w;

    if show_plot && n <= Mplot
        plot3(X_true(1,:), X_true(2,:), X_true(3,:), 'b-',  'LineWidth', 1.2);
        plot3(X_pred(1,:), X_pred(2,:), X_pred(3,:), 'r--', 'LineWidth', 1.2);
    end
end
if show_plot, legend('True','Pred'); end

% ==== 原有 RMSE 汇总 ====
STATS.RMSE_pos_mean  = mean(rmse_pos);
STATS.RMSE_pos_std   = std(rmse_pos);
STATS.RMSE_full_mean = mean(rmse_full);
STATS.RMSE_full_std  = std(rmse_full);
fprintf('\n位置  RMSE: %.4f ± %.4f m\n',  STATS.RMSE_pos_mean,  STATS.RMSE_pos_std);
fprintf('全状态 RMSE: %.4f ± %.4f (量纲化)\n', STATS.RMSE_full_mean, STATS.RMSE_full_std);

% ==== 论文 nRMSE% ====
nrmse_overall.p = 100 * sqrt( safe_div(num_all.p, den_all.p) );
nrmse_overall.v = 100 * sqrt( safe_div(num_all.v, den_all.v) );
nrmse_overall.theta = 100 * sqrt( safe_div(num_all.t, den_all.t) );
nrmse_overall.omega = 100 * sqrt( safe_div(num_all.w, den_all.w) );
vals = [nrmse_overall.p, nrmse_overall.v, nrmse_overall.theta, nrmse_overall.omega];
nrmse_overall.average = mean(vals);

nrmse_traj_mean.p     = mean(traj_pct.p);   nrmse_traj_std.p     = std(traj_pct.p);
nrmse_traj_mean.v     = mean(traj_pct.v);   nrmse_traj_std.v     = std(traj_pct.v);
nrmse_traj_mean.theta = mean(traj_pct.t);   nrmse_traj_std.theta = std(traj_pct.t);
nrmse_traj_mean.omega = mean(traj_pct.w);   nrmse_traj_std.omega = std(traj_pct.w);
avg_vec = [nrmse_traj_mean.p, nrmse_traj_mean.v, nrmse_traj_mean.theta, nrmse_traj_mean.omega];
avg_std = [nrmse_traj_std.p,  nrmse_traj_std.v,  nrmse_traj_std.theta,  nrmse_traj_std.omega ];
nrmse_traj_mean.average = mean(avg_vec);
nrmse_traj_std.average  = mean(avg_std);

STATS.nrmse_pct_overall    = nrmse_overall;
STATS.nrmse_pct_traj.mean  = nrmse_traj_mean;
STATS.nrmse_pct_traj.std   = nrmse_traj_std;
STATS.K_eval = numel(X_list);

fprintf('\n===== nRMSE%% （论文口径）=====\n');
fprintf('Overall（全部样本合并）：\n');
fprintf('  p: %.2f%%,  v: %.2f%%,  Θ: %.2f%%,  ω: %.2f%%,  avg: %.2f%%\n', ...
    nrmse_overall.p, nrmse_overall.v, nrmse_overall.theta, nrmse_overall.omega, nrmse_overall.average);
fprintf('Per-trajectory（均值 ± 标准差）：\n');
fprintf('  p: %.2f ± %.2f %%\n', nrmse_traj_mean.p,     nrmse_traj_std.p);
fprintf('  v: %.2f ± %.2f %%\n', nrmse_traj_mean.v,     nrmse_traj_std.v);
fprintf('  Θ: %.2f ± %.2f %%\n', nrmse_traj_mean.theta, nrmse_traj_std.theta);
fprintf('  ω: %.2f ± %.2f %%\n', nrmse_traj_mean.omega, nrmse_traj_std.omega);
fprintf('  avg: %.2f ± %.2f %%\n', nrmse_traj_mean.average, nrmse_traj_std.average);

%% ===== R² =====


% R²（位置 & 全状态；基于拼接 rollout）
Yt_pos = [];  Yp_pos = [];
Yt_all = [];  Yp_all = [];
for n=1:numel(X_list)
    Xt = X_list{n}; Xp = pred_store{n};
    Yt_pos = [Yt_pos, Xt(1:3,2:end)];
    Yp_pos = [Yp_pos, Xp(1:3,2:end)];
    Yt_all = [Yt_all, Xt(:,2:end)];
    Yp_all = [Yp_all, Xp(:,2:end)];
end
STATS.R2_pos  = r2_score(Yt_pos, Yp_pos);
STATS.R2_full = r2_score(Yt_all, Yp_all);
fprintf('R² (pos):  %.4f\n', STATS.R2_pos);
fprintf('R² (full): %.4f\n', STATS.R2_full);

end % ====== END OF MAIN FUNCTION ======

% ---------- local functions ----------
function z = local_make_z_guard(x, nPhi, strip_const)
    x = x(:);
    if numel(x) ~= 12, error('x 必须为 12×1'); end
    phi = get_basis_x(x); phi = phi(:);
    if strip_const && ~isempty(phi), phi(1) = []; end
    if numel(phi) < nPhi
        phi = [phi; zeros(nPhi-numel(phi),1)];
    elseif numel(phi) > nPhi
        if strip_const && numel(phi) == nPhi+1, phi = phi(2:end);
        else,                                   phi = phi(1:nPhi);
        end
    end
    z = [x; phi];
end

function y = safe_div(a,b)
    if b<=0, y = 0; else, y = a/b; end
end

% R²
function R2 = r2_score(Ytrue, Ypred)
    Ytrue = Ytrue(:); Ypred = Ypred(:);
    SS_res = sum((Ytrue - Ypred).^2);
    SS_tot = sum((Ytrue - mean(Ytrue)).^2) + 1e-12;
    R2 = 1 - SS_res/SS_tot;
end
