function [alpha_1, alpha_2, U_1, U_2, U_3, U_4] = Til_inputU(a, b, c, d, e, f)
   % 求解 U_1
    U_1 = sqrt(a^2 + c^2);
    % 求解 U_2
    U_2 = sqrt(b^2 + d^2);
    % 求解 alpha_1
    alpha_1 = atan2(a, c);
    % 求解 U_3
    U_3 = e;
    % 求解 alpha_2
    alpha_2 = atan2(b, d);
    % 求解 U_4
    U_4 = f;
  
end 