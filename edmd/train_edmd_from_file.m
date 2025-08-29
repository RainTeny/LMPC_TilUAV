%% train_edmd_from_file.m — 从最新数据训练 EDMD（脚本）
clear; clc;

% —— 固定到工程根 —— 
thisFile = mfilename('fullpath');
projRoot = fileparts(thisFile);
cd(projRoot);

% 路径（确保能找到 edmd/ 下的函数）
addpath('C:\Users\Huang_Yuteng\Desktop\MatlabWork\koopmanmpc_quadrotor\DLKoopmanMPC_TilQ\edmd');

%% 1) 载入最新的数据文件（工程根下 data/）
dataDir = fullfile(projRoot,'data');
L = dir(fullfile(dataDir,'dataset_train_*.mat'));
assert(~isempty(L), 'data/ 下没有 dataset_train_*.mat，请先运行 make_dataset.m');
[~,idx] = max([L.datenum]);
dataFile = fullfile(L(idx).folder, L(idx).name);

D = dir(dataFile);
fprintf('使用最新数据文件：%s\n  大小=%.1f MB  时间=%s\n', dataFile, D.bytes/1e6, D.date);

S = load(dataFile, 'X1','X2','U1','t_grid');
X1 = S.X1;  X2 = S.X2;  U1 = S.U1;  t_grid = S.t_grid;

Tminus1  = numel(t_grid)-1;
traj_cnt = size(X1,2) / Tminus1;

fprintf('X1 尺寸：%d × %d（轨迹条数=%d，每条列数=%d）\n', ...
    size(X1,1), size(X1,2), traj_cnt, Tminus1);
fprintf('X2 尺寸：%d × %d\n', size(X2,1), size(X2,2));
fprintf('U1 尺寸：%d × %d\n', size(U1,1), size(U1,2));

assert(size(X1,1)==12 && size(X2,1)==12 && size(U1,1)==6, '维度不符(12/12/6)');
assert(size(X1,2)==size(X2,2) && size(X1,2)==size(U1,2), '样本数不一致');

%% 2) 训练（保持你原来开关；仅新增位置加权）
opts = struct();
opts.use_bias        = true;
opts.use_cross       = true;      % 使用交叉项 (u⊗z)
opts.cross_with_phys = true;      % 传入 X_current；你当前 build_cross 的物理附加仍是注释，不增维
opts.use_norm        = true;      % z-score + 岭回归
opts.reg_scale       = 1e-9;      % 需要更稳可调大：2e-4, 5e-4, 1e-3
opts.strip_const     = true;
opts.block_size      = 6e9;

% === 新增：位置输出加权（更稳更小的 position RMSE）===
opts.pos_weight      = 1.0;         % 给 [x y z] 三行更高权重；不想用就注释掉或设 []

fprintf('\n==> 训练 EDMD 双线性模型 …\n');
fprintf('  use_cross=%d, cross_with_phys=%d, use_norm=%d, strip_const=%d, pos_weight=%s\n', ...
    opts.use_cross, opts.cross_with_phys, opts.use_norm, opts.strip_const, mat2str(opts.pos_weight));
fprintf('  reg_scale=%.1e, block_size=%g\n', opts.reg_scale, opts.block_size);

EDMD = get_EDMD(X1, X2, U1, t_grid, opts);
fprintf('  A: %dx%d   B: %dx%d   E: %dx%d\n', size(EDMD.A,1), size(EDMD.A,2), size(EDMD.B,1), size(EDMD.B,2), size(EDMD.E,1), size(EDMD.E,2));
fprintf('  z 维度 m = %d,  cond(G+R) = %.2e\n', EDMD.nZ, EDMD.condG);


%% 3) 保存模型（同目录 data/，带时间戳）
ts = datestr(now,'yyyymmdd_HHMMSS');
mdlFile = fullfile(dataDir, ['EDMD_model_trained_' ts '.mat']);
save(mdlFile, 'EDMD','-v7.3');
fprintf('√ 模型已保存到：%s\n', mdlFile);

%% 4) 基本自检
assert(size(EDMD.C,1)==12 && size(EDMD.C,2)==size(EDMD.A,1), 'C 尺寸与 A 不一致');
if any(any(EDMD.C(:,1:12) ~= eye(12)))
    warning('C 前 12 列不是单位阵，检查 C 构造。');
end
%% === Export EDMD -> MPC params (robust) ===
% 需要 EDMD 里已有字段：A, B, (可选) b 或 c, (可选) Bj 或 E, t_grid

A = EDMD.A;                          % nz x nz
% ---- Bu 朝向自适配 ----
if size(EDMD.B,1) == size(A,1)
    Bu = EDMD.B;                     % 64x6 已是 (nz x nu)
else
    Bu = EDMD.B.';                   % 6x64 -> 64x6
end
% ---- 常量项 ----
if isfield(EDMD,'c')
    c = EDMD.c;  if isrow(c), c = c.'; end
elseif isfield(EDMD,'b')
    c = EDMD.b;  if isrow(c), c = c.'; end
else
    c = zeros(size(A,1),1);
end
% ---- M（kron(u,z) 的映射）----
if isfield(EDMD,'E')                 % 你日志里就是 64x384
    M = EDMD.E;                      % nz x (nu*nz)
elseif isfield(EDMD,'Bj')
    M = cat(2, EDMD.Bj{:});          % [B1, B2, ... , Bnu]
else
    error('EDMD 里缺少 E 或 Bj 来构造 M');
end

% ---- 采样时间 ----
if isfield(EDMD,'t_grid')
    tg = EDMD.t_grid(:);  Ts = tg(2)-tg(1);
else
    Ts = 0.02; % 没有就先给个默认，按你的实际采样改
end

% ---- 维度（不要再用 m 当两个意思！）----
nz = size(A,1);
nu = size(Bu,2);

% ---- 一致性自检（应接近 0）----
z = randn(nz,1); 
u = randn(nu,1);
z1 = A*z + M*kron(u,z) + Bu*u + c;

% 用 Σ u_j B_j z 重建（若没有 Bj，就把 M 切块）
if isfield(EDMD,'Bj')
    Bj = EDMD.Bj;
else
    Bj = cell(1,nu);
    for j = 1:nu
        cols = (j-1)*nz + (1:nz);
        Bj{j} = M(:, cols);
    end
end
sumterm = zeros(nz,1);
for j=1:nu
    sumterm = sumterm + u(j)*Bj{j}*z;
end
z2 = A*z + sumterm + Bu*u + c;
fprintf('[EDMD->MPC] kron一致性误差 = %.3e\n', norm(z1 - z2));

% ---- 字典中原状态的下标（按你的实际修改）----
dict_info.pos_idx = 1:3;     % x,y,z
dict_info.vel_idx = 4:6;     % dx,dy,dz
dict_info.att_idx = 7:8;     % phi, theta
% 若需要：dict_info.yaw_idx = 9; dict_info.angvel_idx = 10:12;

% ---- 目标保存目录：工程根\mpc ----
thisfile = mfilename('fullpath');   % ...\edmd\tools\train_edmd_from_file.m
root = fileparts(thisfile);
% 向上找直到发现 mpc/ 或到达根目录
while exist(fullfile(root,'mpc'),'dir') ~= 7
    parent = fileparts(root);
    if strcmp(parent, root), break; end
    root = parent;
end
dest_dir = fullfile(root,'mpc');    % 优先根\mpc
if exist(dest_dir,'dir') ~= 7
    mkdir(dest_dir);
end

% 保存
dst_latest = fullfile(dest_dir,'edmd_params.mat');
ts = datestr(now,'yyyymmdd_HHMMSS');
dst_stamp  = fullfile(dest_dir, ['edmd_params_' ts '.mat']);
save(dst_latest, 'A','M','Bu','c','Ts','dict_info','-v7.3');
copyfile(dst_latest, dst_stamp);
fprintf('[EDMD->MPC] 已写入: %s (Ts=%.6f s)\n', dst_latest, Ts);
