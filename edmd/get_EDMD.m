function EDMD = get_EDMD(X1, X2, U1, t_grid, opts)
% GET_EDMD  由 (X1,X2,U1) 训练 Koopman-EDMD 线性模型（两遍流式，省内存）
% 目标保持：z_{k+1} ≈ A z_k + B u_k + E (u_k ⊗ z_k) + b
% —— 最小改动修复：
%   • Bias 行不做标准化（否则被标准成0→G奇异）
%   • 零方差行 std=1 以避免 NaN/Inf
%   • G+R 条件数差时自动加 jitter
%   • 保留 B/E，同时新增 Bj={B_j}
%   • %%% NEW: 训练时对 Z2（输出）加行权重（位置更稳）

% ---------- 1) 默认参数 ----------
if nargin < 5 || isempty(opts), opts = struct(); end
opts = set_default(opts,'use_bias',        true);
opts = set_default(opts,'use_cross',       true);
opts = set_default(opts,'cross_with_phys', true);
opts = set_default(opts,'use_norm',        true);
opts = set_default(opts,'reg_scale',       1e-7);
opts = set_default(opts,'strip_const',     true);
opts = set_default(opts,'block_size',      6e9);
opts = set_default(opts,'params', []);   % 传入动力学参数结构，供 build_cross 物理项使用
opts = set_default(opts,'row_weights', []);   % (可选) 长度 = nZ 的行权重向量；为空则不加权
opts = set_default(opts,'pos_weight',  []);   % (可选) 仅给位置[x,y,z]一个统一权重 w，比如 5

% ---------- 2) 基本检查 ----------
assert(size(X1,1)==12 && size(X2,1)==12 && size(U1,1)==6, '维度不符(12/12/6)');
assert(size(X1,2)==size(X2,2) && size(X1,2)==size(U1,2), '样本数不一致');

N  = size(X1,2);
nU = size(U1,1);
bs = opts.block_size;

% ---------- 3) 字典维度探测 ----------
phi1 = get_basis_x(X1(:,1)); phi1 = phi1(:);
if opts.strip_const && ~isempty(phi1), phi1(1) = []; end
nPhi = numel(phi1);
nZ   = 12 + nPhi;

% 交叉项列数
if opts.use_cross
    z1  = [X1(:,1); phi1];
    if opts.cross_with_phys
        zu1 = build_cross(z1, U1(:,1), X1(:,1));  % 你的 build_cross 的“物理项”当前是注释，不增维
    else
        zu1 = build_cross(z1, U1(:,1), []);
    end
    nZU = numel(zu1);
else
    nZU = 0;
end

% 回归自变量行数
rAug = nZ + nU + nZU + double(opts.use_bias);

%%% NEW: 构造 Z2（输出）行权重（仅作用在目标侧），用 √w 缩放
Wrows = [];  % 空=不加权
if ~isempty(opts.row_weights)
    assert(numel(opts.row_weights)==nZ, 'row_weights 长度应为 nZ=%d', nZ);
    Wrows = sqrt( max(1e-12, opts.row_weights(:)) );
elseif ~isempty(opts.pos_weight)
    Wrows = ones(nZ,1);
    Wrows(1:3) = sqrt( max(1, opts.pos_weight) );  % z 的前 3 行对应 x,y,z
end

% ---------- 4) Pass-1：统计均值/方差 ----------
mu_aug = zeros(rAug,1);  s2_aug = zeros(rAug,1);
mu_Z2  = zeros(nZ,1);    s2_Z2  = zeros(nZ,1);

for s = 1:bs:N
    e  = min(N, s+bs-1);
    nb = e - s + 1;

    Z1_blk = zeros(nZ, nb);
    Z2_blk = zeros(nZ, nb);
    for k = 1:nb
        x1 = X1(:,s+k-1);  x2 = X2(:,s+k-1);
        p1 = get_basis_x(x1); if opts.strip_const && ~isempty(p1), p1(1) = []; end
        p2 = get_basis_x(x2); if opts.strip_const && ~isempty(p2), p2(1) = []; end
        Z1_blk(:,k) = [x1; p1(:)];
        Z2_blk(:,k) = [x2; p2(:)];
    end
    U_blk = U1(:,s:e);

    if opts.use_cross
        ZU_blk = zeros(nZU, nb);
        for k = 1:nb
            if opts.cross_with_phys
                ZU_blk(:,k) = build_cross(Z1_blk(:,k), U_blk(:,k), X1(:,s+k-1));
            else
                ZU_blk(:,k) = build_cross(Z1_blk(:,k), U_blk(:,k), []);
            end
        end
    else
        ZU_blk = [];
    end

    if opts.use_bias
        Z1_aug_blk = [Z1_blk; U_blk; ZU_blk; ones(1,nb)];
    else
        Z1_aug_blk = [Z1_blk; U_blk; ZU_blk];
    end

    mu_aug = mu_aug + sum(Z1_aug_blk, 2);
    s2_aug = s2_aug + sum(Z1_aug_blk.^2, 2);
    mu_Z2  = mu_Z2  + sum(Z2_blk, 2);
    s2_Z2  = s2_Z2  + sum(Z2_blk.^2, 2);
end

mu_aug  = mu_aug / N;
mu_Z2   = mu_Z2  / N;
std_aug = sqrt(max(s2_aug/N - mu_aug.^2, 0)) + 1e-12;
std_Z2  = sqrt(max(s2_Z2 /N - mu_Z2 .^2, 0)) + 1e-12;

% ---- FIX：Bias 行不做标准化（保持恒为1，避免G奇异）----
if opts.use_bias
    idx_bias = rAug;            % 最后一行是 bias
    mu_aug(idx_bias)  = 0;
    std_aug(idx_bias) = 1;
end

% ---------- 5) Pass-2：累加 Ahat/G ----------
Ahat_n = zeros(nZ,  rAug);
G_n    = zeros(rAug, rAug);

for s = 1:bs:N
    e  = min(N, s+bs-1);
    nb = e - s + 1;

    Z1_blk = zeros(nZ, nb);
    Z2_blk = zeros(nZ, nb);
    for k = 1:nb
        x1 = X1(:,s+k-1);  x2 = X2(:,s+k-1);
        p1 = get_basis_x(x1); if opts.strip_const && ~isempty(p1), p1(1) = []; end
        p2 = get_basis_x(x2); if opts.strip_const && ~isempty(p2), p2(1) = []; end
        Z1_blk(:,k) = [x1; p1(:)];
        Z2_blk(:,k) = [x2; p2(:)];
    end
    U_blk = U1(:,s:e);

    if opts.use_cross
        ZU_blk = zeros(nZU, nb);
        for k = 1:nb
            if opts.cross_with_phys
                ZU_blk(:,k) = build_cross(Z1_blk(:,k), U_blk(:,k), X1(:,s+k-1));
            else
                ZU_blk(:,k) = build_cross(Z1_blk(:,k), U_blk(:,k), []);
            end
        end
    else
        ZU_blk = [];
    end

    if opts.use_bias
        Z1_aug_blk = [Z1_blk; U_blk; ZU_blk; ones(1,nb)];
    else
        Z1_aug_blk = [Z1_blk; U_blk; ZU_blk];
    end

    Z1n_blk = (Z1_aug_blk - mu_aug) ./ std_aug;
    Z2n_blk = (Z2_blk     - mu_Z2 ) ./ std_Z2;

    %%% NEW: 输出行加权（只缩放目标侧 Z2）
    if ~isempty(Wrows)
        Z2n_blk = Wrows .* Z2n_blk;  % 每一行乘 √w
    end

    Ahat_n = Ahat_n + Z2n_blk * Z1n_blk.';
    G_n    = G_n    + Z1n_blk * Z1n_blk.';
end

m      = N;
Ahat_n = Ahat_n / m;
G_n    = G_n    / m;

% ---------- 6) 正则化并保证可解 ----------
lambda = opts.reg_scale * trace(G_n)/max(size(G_n,1),1);
G_reg  = G_n + lambda*eye(size(G_n));

% ---- FIX：若 G_reg 条件数过差，自动加 jitter ----
rc = rcond(G_reg);
if ~isfinite(rc) || rc < 1e-12
    jitter = max(1e-8, 1e-6 * mean(diag(G_reg)));
    G_reg  = G_reg + jitter*eye(size(G_reg));
    warning('get_EDMD:illcond', 'G+R 条件数差，已自动加 jitter=%.2e', jitter);
end

K_n    = Ahat_n / G_reg;
condG  = cond(G_reg);

% ---------- 7) 反标准化 ----------
S2 = diag(std_Z2);
S1 = diag(1./std_aug);
K  = S2 * K_n * S1;
c0 = mu_Z2 - S2 * K_n * (mu_aug./std_aug);

% ---------- 8) 拆分 K → A, B, (E), b ----------
off = 0;
EDMD = struct();
EDMD.A = K(:, off+1 : off+nZ);  off = off + nZ;
EDMD.B = K(:, off+1 : off+nU);  off = off + nU;    % 保留 B 字段（线性 U 块）

if opts.use_cross
    EDMD.E = K(:, off+1 : off+nZU);  off = off + nZU;
else
    EDMD.E = [];
end

if opts.use_bias
    EDMD.b = K(:, off+1) + c0;
else
    EDMD.b = c0;
end

% === NEW：把 E 按输入通道重排成 {B_j}（双线性便捷形式） ===
EDMD.Bj = {};
if opts.use_cross && ~isempty(EDMD.E)
    Bj = cell(nU,1);
    for j = 1:nU
        cols = (j-1)*nZ + (1:nZ);
        Bj{j} = EDMD.E(:, cols);           % 每个 m×m
    end
    EDMD.Bj = Bj;
end

% 解码矩阵（前12维直接输出原状态）
C = zeros(12, nZ); C(:,1:12) = eye(12);
EDMD.C = C;

% 元数据 & 标准化参数（供预测复用）
EDMD.opts    = opts;
EDMD.t_grid  = t_grid;
EDMD.lambda  = lambda;
EDMD.condG   = condG;
EDMD.nZ      = nZ;
EDMD.nPhi    = nPhi;
EDMD.nZU     = nZU;
EDMD.K       = K;

EDMD.mu_W  = mu_aug;    % ← W 的均值（含 U/UZ/bias）
EDMD.std_W = std_aug;   % ← W 的标准差（bias 行为 1）
EDMD.mu_Z2 = mu_Z2;
EDMD.std_Z2= std_Z2;

%%% NEW: 记录权重（方便复现实验/日志）
if ~isempty(Wrows)
    EDMD.row_weights = Wrows.^2;   % 保存的是 w，而不是 √w
else
    EDMD.row_weights = [];
end

if any(~isfinite(K(:)))
    warning('EDMD:NaNInf', 'K 中出现 NaN/Inf，请增大正则或检查标准化/数据。');
end
end  % ← 主函数到此结束

% ---------- local helper ----------
function s = set_default(s, field, value)
    if ~isfield(s,field) || isempty(s.(field))
        s.(field) = value;
    end
end

