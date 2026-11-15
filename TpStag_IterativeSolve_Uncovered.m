
clc; close all;

% Define 
vwind = 20 * 1000 / (60^2);    % Wind speed [m/s]
G = 436;                       % Solar irradiance [W/m^2]
Tamb = 17 + 273.15;            % Ambient temp [K]
Tsky = 15 + 273.15;            % Sky temp [K] 
h_wind = 2.8 + 3 * vwind;      % Wind convection coefficient [W/(m^2·K)]
alpha = 0.9;                   % Plate/PV absorptivity
eps_p = 0.9;                   % Plate emissivity
sigma = 5.67e-8;               % Stefan-Boltzmann constant [W/(m^2·K^4)]

% PV Electrical Parameters
eta_ref = 0.25;                % Reference electrical efficiency 
beta_ref = 0.0026;             % Temperature coefficient of efficiency [1/K] (0.26%/K)
T_ref = 25 + 273.15;           % Reference temperature (25 °C in K)

% Temperature-Dependent Electrical Efficiency Function 
eta_el = @(Tp) eta_ref * (1 - beta_ref * (Tp - T_ref));

% Stagnation Residual Function 
% Net Solar Heat Input [G * (alpha - eta_el)] = Thermal Losses [Convection + Radiation]
stagnation_eqn = @(Tp) (G * (alpha - eta_el(Tp))) - ...
                       (h_wind * (Tp - Tamb) + eps_p * sigma * (Tp^4 - Tsky^4));

% Solve Using fzero 
Tp_stag_guess = Tamb + 15;     % Initial guess in Kelvin
options = optimset('Display', 'iter');
Tp_stag_solved = fzero(stagnation_eqn, Tp_stag_guess, options);

eta_stag = eta_el(Tp_stag_solved);


fprintf('\n--- PVT Stagnation Results (with Electrical Generation) ---\n');
fprintf('Stagnation Temperature (T_stag): %.2f K (%.2f °C)\n', Tp_stag_solved, Tp_stag_solved - 273.15);
fprintf('Electrical Efficiency at Stagnation: %.2f%%\n', eta_stag * 100);