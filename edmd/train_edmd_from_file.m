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
opts.reg_scale       = 1e-7;      % 需要更稳可调大：2e-4, 5e-4, 1e-3
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


