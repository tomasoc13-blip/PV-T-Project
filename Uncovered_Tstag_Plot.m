clc; close all;

% --- Define static parameters ---
D_i = 20e-3;                   % Inner diameter [m]
vwind = 20 * 1000 / (60^2);    % Convert km/h to m/s
ro = 998;                      % Density of water [kg/m^3]
mu = 0.890e-3;                 % Dynamic viscosity [Pa·s]
k_water = 0.6;                 % Thermal conductivity of water [W/(m·K)]

params.G = 436;                % Solar irradiance [W/m^2]
params.Tamb = 17 + 273.15;     % Ambient temp [K]
params.Tsky = 15 + 273.15;     % Sky temp [K] 
params.h_wind = 2.8 + 3 * vwind;
params.Tin = 10 + 273.15;      % Inlet temp [K]
params.Cp = 4182;              % Specific heat [J/(kg·K)]
params.Ac = 2;               % Collector area [m^2]
params.alpha = 0.9;            % Absorptivity
params.eps_p = 0.9;            % Emissivity
params.P = 3.14 * D_i;         % Wetted perimeter [m]
params.L = 10;                 % Length [m]

% --- Sweep Flow Rate in L/min (0.001 to 20 L/min) ---
V_lmin_vals = linspace(0.001, 10, 5000); 

Tp_arr = zeros(size(V_lmin_vals));
Q_arr = zeros(size(V_lmin_vals));

options = optimset('Display', 'off'); % Suppress fzero iteration text during sweep

for i = 1:length(V_lmin_vals)
    V_litre_min = V_lmin_vals(i);
    mdot = V_litre_min * ro / 60000; % Convert L/min to kg/sec
    params.mdot = mdot;
    
    % Reynolds and Nusselt calculations
    Re_D = 4 * mdot / (3.14 * D_i * mu);
    Pr = 4182 * mu / k_water;
    
    if Re_D < 2300
        Nu_D = 3.66; % Laminar (constant wall temp assumption)
    elseif Re_D < 3000
        f_3000 = (0.790 * log(3000) - 1.64)^(-2);
        Nu_3000 = ((f_3000 / 8) * (3000 - 1000) * Pr) / ...
                  (1 + 12.7 * sqrt(f_3000 / 8) * (Pr^(2/3) - 1));
        gamma = (Re_D - 2300) / (3000 - 2300);
        Nu_D = (1 - gamma) * 3.66 + gamma * Nu_3000; % Transition interpolation
    else
        f = (0.790 * log(Re_D) - 1.64)^(-2);
        Nu_D = ((f / 8) * (Re_D - 1000) * Pr) / ...
             (1 + 12.7 * sqrt(f / 8) * (Pr^(2/3) - 1)); % Turbulent (Gnielinski)
    end
    params.hfi = k_water * Nu_D / D_i;
    
    % Solve for plate temperature
    Tp_guess = 310; 
    Tp_solved = fzero(@(Tp) bare_residual(Tp, params), Tp_guess, options);
    
    % Store results
    Tp_arr(i) = Tp_solved - 273.15; 
    [~, ~, Q_useful_solved] = bare_residual(Tp_solved, params);
    Q_arr(i) = Q_useful_solved;
end

% --- Calculate exact L/min values for Re = 2300 and Re = 3000 ---
mdot_2300 = (2300 * pi * D_i * mu) / 4;
V_lmin_2300 = mdot_2300 * 60000 / ro;

mdot_3000 = (3000 * pi * D_i * mu) / 4;
V_lmin_3000 = mdot_3000 * 60000 / ro;

% --- Plotting Results ---
figure('Color', 'white', 'Position', [100, 100, 950, 550]);

% --- Shade Transition Region (Re 2300 to 3000) in Light Grey FIRST ---
yl_temp = [params.Tin - 273.15 - 2, max(Tp_arr) + 10];
patch([V_lmin_2300, V_lmin_3000, V_lmin_3000, V_lmin_2300], ...
      [yl_temp(1), yl_temp(1), yl_temp(2), yl_temp(2)], ...
      [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on;

yyaxis left
h_Tp = plot(V_lmin_vals, Tp_arr, 'b-', 'LineWidth', 2);
ylabel('Plate Temperature T_p [°C]', 'FontSize', 12, 'FontWeight', 'bold');
ylim(yl_temp);

yyaxis right
h_Q = plot(V_lmin_vals, Q_arr, '-', 'Color', [0.85, 0.43, 0.1], 'LineWidth', 2);
ylabel('Useful Heat Captured Q_{useful} [W]', 'FontSize', 12, 'FontWeight', 'bold');


xlabel('Flow Rate V (L/min)', 'FontSize', 12, 'FontWeight', 'bold');
title('PV/T T_p and Q_{useful} vs. Flow Rate with Stagnation & Transition Limits', 'FontSize', 13, 'FontWeight', 'bold');

% --- Add reference lines with larger, bold labels ---
h_re2300 = xline(V_lmin_2300, 'k--', 'Re = 2300', 'LineWidth', 1.2, 'LabelVerticalAlignment', 'top', 'FontSize', 11, 'FontWeight', 'bold');
h_re3000 = xline(V_lmin_3000, 'k--', 'Re = 3000', 'LineWidth', 1.2, 'LabelVerticalAlignment', 'top', 'FontSize', 11, 'FontWeight', 'bold');

yyaxis left 
h_stag = yline(29.32, 'b:', 'Stagnation Temp T_{stag} = 29.32 °C', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
hold off;

grid off; 
box off;

% Make axes lines and tick labels bold and prominent
set(gca, 'LineWidth', 1.5, 'FontSize', 11, 'FontWeight', 'bold');

% --- Comprehensive Legend ---
legend([h_Tp, h_Q, h_re2300, h_re3000, h_stag], ...
       {'Plate Temperature (T_p)', 'Useful Heat (Q_{useful})', 'Laminar Limit (Re = 2300)', 'Turbulent Start (Re = 3000)', 'Stagnation Temp (T_{stag})'}, ...
       'Location', 'east', 'FontSize', 10);