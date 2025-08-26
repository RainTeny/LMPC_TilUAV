function zu = build_cross(z, u, X_current, params)
% 控制-状态双线性项：zu = kron(u, z)  [+ 物理项(注释)]
% z: nzx1 ；u: [F1 F2 F3 F4 a1 a2] (弧度)
% X_current: [x y z dx dy dz phi theta psi p q r] (可选)
% params: 可选（B层物理项需要）

    % --- 检查 ---
    if ~isvector(z), error('z 必须是向量'); end
    if ~isvector(u) || numel(u) ~= 6, error('u 必须是长度为6的向量'); end
    z = double(z(:));
    u = double(u(:));

    % --- 主项（双线性必须） ---
    % 注意顺序：kron(u,z) = [u1*z; u2*z; ... u6*z]
    zu = kron(u, z);

    % =========================
    % ====  B层物理项(开)  ====
    % =========================
    if nargin>=3 && ~isempty(X_current) && nargin>=4 && ~isempty(params)
        xcur = double(X_current(:));
        ph=xcur(7); th=xcur(8);
        F1=u(1); F2=u(2); F3=u(3); F4=u(4);
        a1=u(5); a2=u(6);
    
        m = params.mass; g = 9.81;
        Ix=params.J(1,1); Iy=params.J(2,2); Iz=params.J(3,3);
        lm=params.lm; d=params.d;
    
        Fsum_n = (F1+F2+F3+F4)/(m*g);
        Fz_n   = (cos(a1)*F1 + cos(a2)*F3 + F2 + F4)/(m*g);
        Fx_n   = (sin(a1)*F1 + sin(a2)*F3)/(m*g);
        Fz_cth = Fz_n*cos(th);
        Fx_cth = Fx_n*cos(th);
        d13_n  = (F1-F3)/(m*g);
        d24_n  = (F2-F4)/(m*g);
    
        Mx_n = (lm*(cos(a1)*F1 - cos(a2)*F3))/Ix;
        My_n = (lm*(F3 - F4))/Iy;
        Mz_n = (lm*(sin(a1)*F1 + sin(a2)*F3) + d*(cos(a1)*F1 + cos(a2)*F3 - F3 - F4))/Iz;
    
        phys = [Fsum_n; Fz_n; Fx_n; Fz_cth; Fx_cth; d13_n; d24_n; Mx_n; My_n; Mz_n];
        zu = [zu; phys(:)];
    end
end

