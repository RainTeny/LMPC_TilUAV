%% eval_edmd_from_file.m — 自动加载最新模型并评估（含论文 nRMSE% 表）
clear; clc;

% ---------- 1) 固定到脚本所在目录（工程根/edmd 目录） ----------
thisFile = mfilename('fullpath');
here     = fileparts(thisFile);
cd(here);

addpath(here);
if exist(fullfile(here,'data'),'dir'), addpath(fullfile(here,'data')); end

% ---------- 2) 找到 data/ 下最新的模型文件 ----------
dataDir = fullfile(here,'data');
assert(exist(dataDir,'dir')==7, '未找到数据目录：%s', dataDir);

Lm = dir(fullfile(dataDir,'EDMD_model_trained_*.mat'));
if isempty(Lm)
    Lm = dir(fullfile(dataDir,'EDMD_model_trained.mat'));
end
assert(~isempty(Lm), 'data/ 下没有 EDMD 模型文件，请先训练。');
[~,idxm] = max([Lm.datenum]);
modelFile = fullfile(Lm(idxm).folder, Lm(idxm).name);
Dm = dir(modelFile);
fprintf('使用最新模型：%s\n  大小=%.1f MB  时间=%s\n', modelFile, Dm.bytes/1e6, Dm.date);

S = load(modelFile,'EDMD'); EDMD = S.EDMD;

% ---------- 3) 评估参数 ----------
dt        = 0.01;
t_span    = 1;
X0        = zeros(12,1);
show_plot = true;
seed      = [];

% ---------- 4) 执行评估 ----------
STATS = eval_EDMD(X0, dt, t_span, EDMD, show_plot, seed);

% ---------- 5) 按论文表格打印 ----------
P = STATS.nrmse_pct_traj.mean;
Sg= STATS.nrmse_pct_traj.std;
O = STATS.nrmse_pct_overall;


