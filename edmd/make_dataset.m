%% make_dataset.m — 采样并保存训练数据（脚本）
clear; clc;

% —— 将工作目录固定到本脚本所在目录（工程根）——
thisFile = mfilename('fullpath');
projRoot = fileparts(thisFile);
cd(projRoot);

%% 采样参数
dt        = 0.001;
t_grid    = 0:dt:0.1;            % ⇒ T = 201, Tminus1 = 200
need_good = 2000;                 % 目标有效轨迹条数
batch_try = 1000;
show_plot = false;
flag_mode = 'train';

Tminus1 = numel(t_grid)-1;       % 每条轨迹的列数 = T-1
X0      = zeros(12,1);

X1 = []; X2 = []; U1 = [];
good_cnt = 0; batch_id = 0;

fprintf('==> 开始采样 …\n');
while good_cnt < need_good
    batch_id = batch_id + 1;

    [~,~,Xi1,Xi2,Ui1,~] = get_rnd_trajectories( ...
        X0, batch_try, t_grid, show_plot, flag_mode);

    n_new    = size(Xi1,2) / Tminus1;   % 本批有效轨迹条数
    good_cnt = good_cnt + n_new;

    X1 = [X1 Xi1];   %#ok<AGROW>
    X2 = [X2 Xi2];
    U1 = [U1 Ui1];

    fprintf('  批 %d：新增 %d 条，累计 %d / %d\n', ...
        batch_id, n_new, good_cnt, need_good);
end
fprintf('采样完成，累计有效轨迹 = %d 条（可能超额）\n', good_cnt);

%% —— 强制裁剪到 need_good 条（防止超额）——
assert(mod(size(X1,2),Tminus1)==0 && mod(size(X2,2),Tminus1)==0 && mod(size(U1,2),Tminus1)==0, ...
    'X1/X2/U1 的列数不是 T-1 的整数倍');

traj_cnt = size(X1,2)/Tminus1;
if traj_cnt > need_good
    idx_traj = 1:need_good;              % 或 randperm(traj_cnt,need_good)
    idx_cols = [];
    for k = idx_traj
        idx_cols = [idx_cols (k-1)*Tminus1 + (1:Tminus1)]; %#ok<AGROW>
    end
    X1 = X1(:,idx_cols);
    X2 = X2(:,idx_cols);
    U1 = U1(:,idx_cols);
    traj_cnt = need_good;
end
fprintf('最终写盘轨迹条数 = %d（每条列数 = %d）\n', traj_cnt, Tminus1);

%% 保存到工程根下 data/
dataDir = fullfile(projRoot,'data');
if ~exist(dataDir,'dir'), mkdir(dataDir); end
ts = datestr(now,'yyyymmdd_HHMMSS');
outFile = fullfile(dataDir, ['dataset_train_' ts '.mat']);
save(outFile, 'X1','X2','U1','t_grid','-v7.3');

info = dir(outFile);
fprintf('√ 训练数据已保存到：%s\n  大小=%.1f MB  时间=%s\n', ...
    outFile, info.bytes/1e6, info.date);

