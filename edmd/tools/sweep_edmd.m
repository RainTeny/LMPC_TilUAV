%% sweep_edmd.m — 强化版横扫：更大参数空间 & 自动输出最佳参数
clear; clc;

%% === 0) 工程与路径（鲁棒 addpath） ===
thisFile = mfilename('fullpath');
if isempty(thisFile), projRoot = pwd; else, projRoot = fileparts(thisFile); end
cd(projRoot);

candPaths = { fullfile(projRoot,'edmd'), fullfile(projRoot,'edmd','tools'), fullfile(projRoot,'tools'), fullfile(projRoot,'src') };
for i=1:numel(candPaths)
    if exist(candPaths{i},'dir'), addpath(candPaths{i}); end
end
needFns = {'get_EDMD','get_basis_x','build_cross'};
for i=1:numel(needFns)
    if exist(needFns{i},'file')~=2
        error('找不到函数 %s，请确认路径设置正确。', needFns{i});
    end
end

%% === 1) 数据定位（更宽松 + 递归搜索） ===
dataDirOverride = '';                % 可填指定文件或目录；留空则自动搜索
pattern = 'dataset_train*.mat';

dataFile = '';
if ~isempty(dataDirOverride)
    if exist(dataDirOverride,'file')==2
        dataFile = dataDirOverride;
    elseif exist(dataDirOverride,'dir')==7
        L = dir(fullfile(dataDirOverride, pattern));
        if ~isempty(L), [~,idx]=max([L.datenum]); dataFile = fullfile(L(idx).folder,L(idx).name); end
    else
        warning('dataDirOverride 不是文件也不是目录：%s', dataDirOverride);
    end
end
if isempty(dataFile)
    dataDir1 = fullfile(projRoot,'data');
    L = [];
    if exist(dataDir1,'dir'), L = dir(fullfile(dataDir1, pattern)); end
    if isempty(L)
        fprintf('在 %s/data 未找到 %s，递归搜索全部子目录...\n', projRoot, pattern);
        Lall = [];
        allPaths = regexp(genpath(projRoot), pathsep, 'split');
        for i=1:numel(allPaths)
            p = allPaths{i}; if isempty(p) || ~isfolder(p), continue; end
            Li = dir(fullfile(p, pattern));
            if ~isempty(Li), Lall=[Lall; Li(:)]; end %#ok<AGROW>
        end
        L = Lall;
    end
    assert(~isempty(L), '未找到任何匹配 %s 的数据集，请先生成数据。', pattern);
    [~,idx] = max([L.datenum]);
    dataFile = fullfile(L(idx).folder, L(idx).name);
end

D = dir(dataFile);
fprintf('使用数据文件：%s\n  大小=%.1f MB  时间=%s\n', dataFile, D.bytes/1e6, D.date);
S = load(dataFile, 'X1','X2','U1','t_grid');
X1 = S.X1;  X2 = S.X2;  U1 = S.U1;  t_grid = S.t_grid;

T      = numel(t_grid);
Tminus1= T-1;
trajN  = size(X1,2) / Tminus1;
assert(mod(size(X1,2),Tminus1)==0, '列数与每轨迹步数不整除');
fprintf('数据：轨迹=%d，单条步数=%d，Δt≈%.4fs\n', trajN, Tminus1, t_grid(2)-t_grid(1));

%% === 2) 训练/验证划分（按轨迹分） ===
val_ratio   = 0.2;
val_trajN   = max(1, round(trajN*val_ratio));
train_trajN = trajN - val_trajN;
train_ids   = 1:train_trajN;
val_ids     = (train_trajN+1):trajN;
col_range   = @(i_traj) ((i_traj-1)*Tminus1 + 1) : (i_traj*Tminus1);
train_cols  = cell2mat(arrayfun(@(i) col_range(i), train_ids, 'UniformOutput', false));
val_cols    = cell2mat(arrayfun(@(i) col_range(i), val_ids,   'UniformOutput', false));

%% === 3) 横扫参数空间（可改） ===
reg_list         = logspace(-6, -1, 12);   % 从 1e-6 到 1e-1（更大范围）
strip_options    = [true, false];          % 先全扫，看看是否“去常数”更好
cphys_options    = [true];                  % 你现在 build_cross 的物理项注释了；如要对比可改为 [true false]

fprintf('横扫 reg_scale = %s\n', mat2str(reg_list,3));
fprintf('横扫 strip_const ∈ {%s}\n', join(string(strip_options), ", "));
fprintf('横扫 cross_with_phys ∈ {%s}\n', join(string(cphys_options), ", "));

% 评估设置
horizons  = [10, 25, 50, 100];  % 0.1 / 0.25 / 0.5 / 1.0 s (Δt=0.01s)
eval_dims = [1 2 3 7 8 9];      % x y z phi theta psi
use_nrmse = true;

%% === 4) 容器 & 迭代 ===
comb = {};
for a = 1:numel(reg_list)
  for b = 1:numel(strip_options)
    for c = 1:numel(cphys_options)
      comb{end+1,1} = [a b c]; %#ok<AGROW>
    end
  end
end
M = numel(comb);
metrics = repmat(struct('reg',[],'strip',[],'cphys',[],'cond',[], ...
                        'nrmse_avg',[],'nrmse_byH',[],'file',''), M, 1);

best_idx = 1; best_score = inf;

%% === 5) 逐组合训练 & 验证 ===
dataDirOut = fileparts(dataFile);
for i = 1:M
    idx = comb{i};
    reg   = reg_list(idx(1));
    strip = strip_options(idx(2));
    cphys = cphys_options(idx(3));

    fprintf('\n==> 训练（reg_scale=%.3e, strip_const=%d, cross_with_phys=%d）\n', reg, strip, cphys);

    opts = struct();
    opts.use_bias        = true;
    opts.use_cross       = true;
    opts.cross_with_phys = cphys;
    opts.use_norm        = true;
    opts.reg_scale       = reg;
    opts.strip_const     = strip;
    opts.block_size      = 1e4;

    EDMD_i = get_EDMD(X1(:,train_cols), X2(:,train_cols), U1(:,train_cols), t_grid, opts);

    [nrmse_byH, nrmse_avg] = eval_on_val(EDMD_i, X1, U1, t_grid, val_ids, horizons, eval_dims, use_nrmse);

    metrics(i).reg       = reg;
    metrics(i).strip     = strip;
    metrics(i).cphys     = cphys;
    metrics(i).cond      = EDMD_i.condG;
    metrics(i).nrmse_avg = nrmse_avg;
    metrics(i).nrmse_byH = nrmse_byH;

    ts = datestr(now,'yyyymmdd_HHMMSS');
    mdlFile = fullfile(dataDirOut, sprintf('EDMD_model_reg%.0e_strip%d_cphys%d_%s.mat', reg, strip, cphys, ts));
    save(mdlFile, 'EDMD_i','-v7.3');
    metrics(i).file = mdlFile;

    fprintf('  cond(G+R)=%.2e  |  nRMSE(10/25/50/100)=[%.3f %.3f %.3f %.3f]  |  avg=%.3f\n', ...
        EDMD_i.condG, nrmse_byH(1), nrmse_byH(2), nrmse_byH(3), nrmse_byH(4), nrmse_avg);

    if isfinite(nrmse_avg) && (nrmse_avg < best_score)
        best_score = nrmse_avg; best_idx = i;
    end
end

%% === 6) 结果汇总、CSV、最佳模型 ===
fprintf('\n===== 横扫结果 =====\n');
Tsum = table('Size',[M 8], ...
    'VariableTypes',{'double','logical','logical','double','double','double','double','string'}, ...
    'VariableNames',{'reg_scale','strip_const','cross_with_phys','cond','nRMSE_H10','nRMSE_H25','nRMSE_H50','model_file'});
for i=1:M
    Tsum.reg_scale(i)     = metrics(i).reg;
    Tsum.strip_const(i)   = metrics(i).strip;
    Tsum.cross_with_phys(i)= metrics(i).cphys;
    Tsum.cond(i)          = metrics(i).cond;
    v = metrics(i).nrmse_byH;
    Tsum.nRMSE_H10(i)     = v(1);
    if numel(v)>=2, Tsum.nRMSE_H25(i)=v(2); end
    if numel(v)>=3, Tsum.nRMSE_H50(i)=v(3); end
    Tsum.model_file(i)    = string(metrics(i).file);
    fprintf('reg=%.3e | strip=%d | cphys=%d | cond=%.2e | nRMSE=[%.3f %.3f %.3f ...] | avg=%.3f\n', ...
        metrics(i).reg, metrics(i).strip, metrics(i).cphys, metrics(i).cond, v(1), v(min(2,end)), v(min(3,end)), metrics(i).nrmse_avg);
end

% 选最优（平均 nRMSE 最小）
best = metrics(best_idx);
fprintf('\n★ 最优参数：reg_scale=%.3e, strip_const=%d, cross_with_phys=%d\n', best.reg, best.strip, best.cphys);
fprintf('  cond(G+R)=%.2e | nRMSE(10/25/50/100)=[%.3f %.3f %.3f %.3f] | avg=%.3f\n', ...
    best.cond, best.nrmse_byH(1), best.nrmse_byH(2), best.nrmse_byH(3), best.nrmse_byH(4), best.nrmse_avg);

bestOut = fullfile(dataDirOut,'EDMD_model_best.mat');
copyfile(best.file, bestOut);
fprintf('√ 已复制最优模型到：%s\n', bestOut);

% 保存 CSV
csvName = fullfile(dataDirOut, sprintf('EDMD_sweep_results_%s.csv', datestr(now,'yyyymmdd_HHMMSS')));
writetable(Tsum, csvName);
fprintf('√ 结果表已保存：%s\n', csvName);

%% ========= 内部辅助函数 =========
function [nrmse_byH, nrmse_avg] = eval_on_val(EDMD, X1, U1, t_grid, val_ids, horizons, eval_dims, use_nrmse)
    Tminus1 = numel(t_grid)-1;
    m       = EDMD.nZ;       
    Hs      = horizons(:).';
    nH      = numel(Hs);
    err_acc = zeros(1,nH);
    cnt_acc = zeros(1,nH);

    for id = val_ids
        cols = ((id-1)*Tminus1 + 1) : (id*Tminus1);
        Xseq = X1(:, cols);
        Useq = U1(:, cols);

        z0 = [ Xseq(:,1);  phi_x_with_strip(Xseq(:,1), EDMD.opts.strip_const) ];

        for h_idx = 1:nH
            H = Hs(h_idx);
            H = min(H, size(Xseq,2)-1);

            z_pred = rollout_bilinear(EDMD, z0, Useq(:,1:H));
            x_hat  = z_pred(1:12);
            x_ref  = Xseq(:,H+1);

            e = x_hat(eval_dims) - x_ref(eval_dims);
            if use_nrmse
                denom = max(1e-8, rms(Xseq(eval_dims,1:H), 2));
                e = e ./ denom;
            end
            err_acc(h_idx) = err_acc(h_idx) + norm(e);
            cnt_acc(h_idx) = cnt_acc(h_idx) + 1;
        end
    end
    nrmse_byH = err_acc ./ max(1,cnt_acc);
    nrmse_avg = mean(nrmse_byH);

    function z_next = rollout_bilinear(EDMD, z0, Useq)
        z = z0;
        for k = 1:size(Useq,2)
            u = Useq(:,k);
            uz_term = zeros(EDMD.nZ,1);
            if ~isempty(EDMD.Bj)
                for j = 1:numel(EDMD.Bj)
                    uz_term = uz_term + u(j) * (EDMD.Bj{j} * z);
                end
            elseif isfield(EDMD,'E') && ~isempty(EDMD.E)
                uz_term = EDMD.E * kron(u, z);
            end
            z = EDMD.A * z + EDMD.B * u + uz_term + EDMD.b;
        end
        z_next = z;
    end

    function p = phi_x_with_strip(x, strip_const_flag)
        p = get_basis_x(x);
        if strip_const_flag && ~isempty(p), p(1)=[]; end
    end
end

