function phi = get_basis_x(X)
% φ_x(x): 只含状态，不含控制
% X: [x y z dx dy dz phi theta psi p q r]  (弧度)

    % --- 检查 ---
    if ~isvector(X) || numel(X) ~= 12
        error('get_basis_x: X 必须是 12x1 或 1x12 的向量');
    end
    x = double(X(:));

    % 解包
    x1=x(1); x2=x(2); x3=x(3);
    dx=x(4); dy=x(5); dz=x(6);
    ph=x(7); th=x(8); ps=x(9);
    p =x(10); q =x(11); r =x(12);

    % 角度包裹
    ph = atan2(sin(ph),cos(ph));
    th = atan2(sin(th),cos(th));
    ps = atan2(sin(ps),cos(ps));

    % ---------- 基础（含位置） ----------
    base = [
        1;                  % 常数
        x1; x2; x3;         % 位置（MPC跟踪解码）
        sin(ph);  cos(ph);  % 姿态周期性
        sin(th);  cos(th);
        sin(ps);  cos(ps);
        dx; dy; dz;         % 线速度（阻尼线性）
        p;  q;  r;          % 角速度（陀螺项）
        dx.^2; dy.^2; dz.^2;% 轻量二次（抗噪）
        p.^2;  q.^2;  r.^2;
        dx.*p;  dy.*q;  dz.*r;
        p.*q;  q.*r;  r.*p
    ];

    % ---------- A层增强：速度×姿态（原来那组，保留） ----------
    v_trig = [
        dx.*cos(th);
        dx.*sin(th);
        dy.*cos(th);
        dy.*sin(th);
        dz.*cos(ph).*cos(th);
        dz.*sin(ph).*cos(th)
    ];

    % ---------- A层增强：角速×姿态（原来那组，保留） ----------
    w_trig = [
        p.*cos(ph);
        p.*sin(ph);
        q.*cos(th);
        q.*sin(th);
        r.*cos(ps);
        r.*sin(ps)
    ];

    % ---------- 新增1：与 yaw 相关的横向耦合（轻量） ----------
    v_yaw = [
        dx.*cos(ps);
        dx.*sin(ps);
        dy.*cos(ps);
        dy.*sin(ps)
    ];

    % ---------- 新增2：速度-角速度的交叉轴耦合（补齐） ----------
    v_w_cross = [
        dx.*q;  dx.*r;
        dy.*p;  dy.*r;
        dz.*p;  dz.*q
    ];

    % ---------- 方向余弦第三列（平动受力方向） ----------
    % R*e3 = [sinθ; -sinφ cosθ; cosφ cosθ]
    dir_cos = [
        sin(th);
        -sin(ph).*cos(th);
        cos(ph).*cos(th)
    ];

    % 汇总（顺序无所谓，但保持稳定）
    phi = [base; v_trig; w_trig; v_yaw; v_w_cross; dir_cos];
    phi = phi(:);
end

