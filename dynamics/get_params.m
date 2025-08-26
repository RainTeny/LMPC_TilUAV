function params = get_params()
% Output
% parameters for the tilt-rotor quadrotor
% Compatible with Euler-angle dynamics model

%% Physical parameters

% Updated values (from your Euler-based tilt rotor model)
m = 0.48;  % mass in kg
J = diag([5.5e-3, 5.5e-3, 9.0e-3]);  % inertia tensor (kg·m^2)
g = 9.81;  % gravity (m/s^2)

lm = 0.18;       % distance from center of mass to rotor (m)
d = 1.574e-7;    % drag coefficient or motor asymmetry

% Damping coefficients (air resistance)
kx = 0.25;
ky = 0.25;
kz = 0.25;

%% Assign to params struct
params.mass = m;
params.J = J;
params.g = g;

params.lm = lm;
params.d = d;

params.kx = kx;
params.ky = ky;
params.kz = kz;
end


