function [X1f, X2f, U1f] = preprocess_dataset(X1, X2, U1, t_grid, opts)
% 轻微平滑 X/U 序列：butterworth 低通 + 零相位 filtfilt
% —— 逐轨迹处理（避免跨轨迹串扰）
% —— 欧拉角先 unwrap，再滤波，最后 wrap 回 [-pi,pi]
%
% 输入
%   X1: 12×N  （时刻 k）
%   X2: 12×N  （时刻 k+1）
%   U1:  6×N  （与 X1 对齐）
%   t_grid: 1×T  （每条轨迹的时间栅格；T-1 = 每条的列数）
%   opts:
%     .method   = 'butter' | 'movavg' | 'sgolay'   （默认 'butter'）
%     .fc_hz    = 10       % 截止频率[Hz]，默认 10（dt=0.01→fs=100Hz，Wn≈0.2）
%     .order    = 2        % 巴特沃斯阶数，默认 2
%     .win      = 7        % 移动平均窗口（奇数）
%     .sgolay_k = 2        % Savitzky-Golay 多项式阶数
%     .sgolay_f = 7        % Savitzky-Golay 帧长（奇数）
%
% 输出
%   X1f, X2f, U1f : 与输入同尺寸的滤波后数据

if nargin<5, opts=struct(); end
if ~isfield(opts,'method'), opts.method='butter'; end
if ~isfield(opts,'fc_hz'),  opts.fc_hz=10;      end
if ~isfield(opts,'order'),  opts.order=2;       end
if ~isfield(opts,'win'),    opts.win=7;         end
if ~isfield(opts,'sgolay_k'), opts.sgolay_k=2;  end
if ~isfield(opts,'sgolay_f'), opts.sgolay_f=7;  end

dt  = t_grid(2)-t_grid(1);
fs  = 1/dt;
T   = numel(t_grid);
Tm1 = T-1;

% 预分配
X1f = zeros(size(X1)); X2f = zeros(size(X2)); U1f = zeros(size(U1));

% 滤波器设计（如用 butter）
if strcmpi(opts.method,'butter')
    Wn = min(0.99, opts.fc_hz/(fs/2));   % 规范化截止频率
    [b,a] = butter(opts.order, Wn, 'low');
end

% 每条轨迹单独处理
trajN = size(X1,2)/Tm1;
assert(mod(size(X1,2),Tm1)==0, '列数与每轨迹步数不整除');

for i = 1:trajN
    idx = (i-1)*Tm1 + (1:Tm1);     % 本轨迹的列索引

    % 还原“整条状态序列”：X_full = [X1(:,k0), X2(:,k0:kend)]
    X_full = [ X1(:,idx(1)), X2(:,idx) ];   % 12×T
    U_full = [ U1(:,idx), U1(:,idx(end)) ]; % 6×T（末端补1列，便于同长滤）

    % 角度 unwrap（φ,θ,ψ = 7,8,9）
    ang_idx = [7 8 9];
    Xang = X_full(ang_idx,:);
    Xang_u = unwrap(Xang, [], 2);  % 沿时间维 unwrap
    X_full(ang_idx,:) = Xang_u;

    % —— 逐通道滤波 —— 
    switch lower(opts.method)
        case 'butter'
            Xf_full = filtfilt(b,a, X_full.').';   % 零相位
            Uf_full = filtfilt(b,a, U_full.').';
        case 'movavg'
            k = opts.win;
            Xf_full = movmean(X_full, k, 2, 'Endpoints','shrink');
            Uf_full = movmean(U_full, k, 2, 'Endpoints','shrink');
        case 'sgolay'
            k = opts.sgolay_k; f = opts.sgolay_f;
            Xf_full = sgolayfilt(X_full, k, f, [], 2);
            Uf_full = sgolayfilt(U_full, k, f, [], 2);
        otherwise
            error('未知 method=%s', opts.method);
    end

    % 角度 wrap 回 [-pi,pi]
    Xf_full(ang_idx,:) = wrapToPi(Xf_full(ang_idx,:));

    % —— 回填到 X1/X2/U1 —— 
    X1f(:,idx) = Xf_full(:,1:end-1);
    X2f(:,idx) = Xf_full(:,2:end);
    U1f(:,idx) = Uf_full(:,1:end-1);
end

% 输入侧的基本合法性保护（仅限训练数据用，不改变真实系统）
% 推力非负
thrust_idx = 1:4;
U1f(thrust_idx,:) = max(U1f(thrust_idx,:), 0);

% 倾转角软限（例如 ±25°），按你当前训练范围改
tilt_idx = 5:6;
U1f(tilt_idx,:) = max(min(U1f(tilt_idx,:), deg2rad(25)), -deg2rad(25));
end
