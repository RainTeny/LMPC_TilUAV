
function [Phi_Phi,Phi_F,Phi_R,A_e,B_e,C_e,Phi_D,F,Phi] = ...
         mpcgain_mimo(A,B,C,Nc,Np,D,Qy)

H = C;
[m1,~] = size(H);  [n_x,nu] = size(B);

A_e = blkdiag(A, eye(m1));
A_e(n_x+1:end,1:n_x) = H*A;
B_e = [B; H*B];
D_e = [D; H*D];
C_e = [zeros(m1, n_x) , eye(m1)];

ny = m1;  nx_aug = n_x + m1;
F = zeros(Np*ny, nx_aug);
Phi = zeros(Np*ny, Nc*nu);
Phi_d = zeros(Np*ny, Nc*size(D_e,2));
for j = 1:Np
    F((j-1)*ny+1:j*ny,:) = C_e*(A_e^j);
    for i = 1:Nc
        if j>=i
            idxY = (j-1)*ny+1 : j*ny;
            idxU = (i-1)*nu+1 : i*nu;
            Phi(idxY, idxU) = C_e*(A_e^(j-i))*B_e;
            Phi_d(idxY,(i-1)*size(D_e,2)+1:i*size(D_e,2)) = ...
                      C_e*(A_e^(j-i))*D_e;
        end
    end
end

Phi_Phi = Phi'*Qy*Phi;
Phi_F   = Phi'*Qy*F;
Phi_D   = Phi'*Qy*Phi_d;
Phi_R   = Phi_F(:, end-ny+1:end);
end
