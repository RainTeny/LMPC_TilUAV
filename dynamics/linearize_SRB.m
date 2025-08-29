function [A_lin, B_lin] = linearize_SRB()
    % Linearize the custom Euler-angle tilt rotor dynamics

    % Symbolic state and input
    x = sym('x', [12, 1]);  % 12 states
    u = sym('u', [6, 1]);   % 6 control inputs
    t = 0;

    % Get system parameters
    params = get_params();

    % Define symbolic dynamics
    f = dynamics_SRB(t, x, u, params);

    % Linearize around hover condition
    x0 = zeros(12, 1);  % zero state
    u0 = zeros(6, 1);   % zero control input: F1-F4=0, alpha1/2=0

    % Jacobians
    A_lin = double(subs(jacobian(f, x), [x; u], [x0; u0]));
    B_lin = double(subs(jacobian(f, u), [x; u], [x0; u0]));
end
