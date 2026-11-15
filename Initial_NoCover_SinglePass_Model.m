clc, close all

% Define inputs
D_i = 20*(10^-3); % Convert mm to m, inner diameter
vwind = 20*(1000)/(60^2); % convert kmh to m/s
ro = 998;
V_litre_min = 0.001; % Design parameter
mdot = V_litre_min*ro/60000; % kg/sec
mu = 0.890*10^-3;
k_water = 0.6;

params.G = 800; params.Tamb = 18 + 273.15; params.Tsky = 10 + 273.15; 
params.h_wind = 2.8 + 3*vwind;
params.Tin = 10 +273.15;
params.mdot = mdot;
params.Cp = 4182;
params.Ac = 2;
params.alpha = 0.9; % Kirchoff's Law
params.eps_p = 0.9;
params.P = 3.14*(D_i); % Wetted perimeter
params.L = 10; % Length

Re_D = 4*mdot/(3.14*D_i*mu);
Pr = 4182*mu/k_water;


if Re_D < 2300
    % Laminar regime for constant wall temperature
    Nu_D = 3.66;
else
    % Transitional / Turbulent regime: Gnielinski correlation
f = (0.790 * log(Re_D) - 1.64)^-2; % Darcy Friction Factor, Flionenko/peutkhov eqn.

    % Compute Gnielinski Nu
    numerator = (f / 8) * (Re_D - 1000) * Pr;
    denominator = 1 + 12.7 * sqrt(f / 8) * (Pr^(2/3) - 1);
    Nu_D = numerator / denominator;
end





params.hfi = k_water*Nu_D/(D_i);


% Initial guess for Tp in Kelvin
Tp_guess = 310; 

% Run fzero
options = optimset('Display', 'iter'); % Shows the solver converging
Tp_solved = fzero(@(Tp) bare_residual(Tp, params), Tp_guess);

fprintf('Solved Plate Temperature: %.2f K (%.2f °C)\n', Tp_solved, Tp_solved - 273.15);

[~, Tout_solved, Q_useful_solved] = bare_residual(Tp_solved, params);

fprintf('Outlet Temp (Tout): %.2f °C\n', Tout_solved - 273.15);
fprintf('Useful Heat (Q_useful): %.2f W\n', Q_useful_solved);