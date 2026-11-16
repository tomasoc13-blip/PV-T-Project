clc, close all, clear

% Define inputs
D_i = 20 * (10^-3); % Convert mm to m, inner diameter
vwind = 20 * (1000) / (60^2); % convert kmh to m/s
ro = 998;
V_litre_min = 2; % Design parameter (L/min)
mdot = V_litre_min * ro / 60000; % kg/sec
mu = 0.890 * 10^-3;
k_water = 0.6;

params.G = 436; 
params.Tamb = 15 + 273.15; 
params.Tsky = 10 + 273.15; 
params.h_wind = 2.8 + 3 * vwind;
params.mdot = mdot;
params.Cp = 4182;
params.Ac = 2;
params.alpha = 0.9; 
params.eps_p = 0.9;
params.P = 3.14 * (D_i); % Wetted perimeter
params.L = 10; % Length

% --- Tank & Euler Time-Stepping Setup ---
V_tank_L = 200; % Tank volume (Litres)
m_tank = (V_tank_L / 1000) * ro; % Tank water mass (kg)
T_tank = 10 + 273.15; % Initial tank temperature (K)
T_target = 45 + 273.15; % Target temperature (45 °C in K)

dt = 60; % Time step: 1 hour (seconds)
time_hrs = 0;
max_hrs = 200; % Safety limit

% Initial guess for Tp in Kelvin
Tp_guess = 310; 
options = optimset('Display', 'off'); % Suppress iteration spam in loop

% Arrays for logging/plotting
hist_time = [0];
hist_Ttank = [T_tank - 273.15];
hist_Quseful = [];
hist_Qelec = [];

fprintf('\n--- GLAZED PVT TANK RECIRCULATION SIMULATION ---\n');
fprintf('Time (hr) | T_tank (°C) | T_plate (°C) | Q_useful (W) | Q_elec (W) | Total Power (W)\n');
fprintf('------------------------------------------------------------------------------------\n');

while T_tank < T_target && time_hrs < max_hrs
    % 1. Set current tank temperature as collector inlet temperature
    params.Tin = T_tank;
    
    % 2. Internal tube flow calculations (Gnielinski / Laminar) based on current Tin properties
    Re_D = 4 * mdot / (3.14 * D_i * mu);
    Pr = 4182 * mu / k_water;
    if Re_D < 2300
        Nu_D = 3.66;
    else
        f = (0.790 * log(Re_D) - 1.64)^-2; 
        numerator = (f / 8) * (Re_D - 1000) * Pr;
        denominator = 1 + 12.7 * sqrt(f / 8) * (Pr^(2/3) - 1);
        Nu_D = numerator / denominator;
    end
    params.hfi = k_water * Nu_D / (D_i);
    
    % 3. Solve steady-state glazed residual for this hour
    Tp_solved = fzero(@(Tp) glazed_residual(Tp, params), Tp_guess, options);
    [~, Tout_solved, Q_useful_solved, h_gap_solved, Q_elec_solved] = glazed_residual(Tp_solved, params);
    
    % Warm-start guess for next step
    Tp_guess = Tp_solved;
    P_pump = 8;
    % Thermal diode guard: prevent reverse heat loss if tank exceeds stagnation plate temp
    Q_thermal_net = max(0, Q_useful_solved);
    Q_total = Q_thermal_net + (Q_elec_solved - P_pump);
    
    % Log current step
    fprintf('%9.1f | %11.2f | %12.2f | %12.2f | %10.2f | %15.2f\n', ...
        time_hrs, T_tank - 273.15, Tp_solved - 273.15, Q_thermal_net, Q_elec_solved, Q_total);
    
    hist_Quseful(end+1) = Q_thermal_net;
    hist_Qelec(end+1) = Q_elec_solved;
    
    % 4. Forward Euler step for fully mixed tank (heated by thermal + electrical)
    dT = (Q_total * dt) / (m_tank * params.Cp);
    T_tank = T_tank + dT;
    time_hrs = time_hrs + (dt / 3600);
    
    hist_time(end+1) = time_hrs;
    hist_Ttank(end+1) = T_tank - 273.15;
end

fprintf('------------------------------------------------------------------------------------\n');
fprintf('Simulation complete! Target 45 °C reached at %.1f hours.\n', time_hrs);

% Plot Tank Temperature vs Time
figure;
plot(hist_time, hist_Ttank, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
yline(45, 'r--', 'Target 45 °C', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Tank Temperature (°C)');
title('200 L Mixed Tank Heating (Glazed Recirculating PVT + Q_{elec})');
grid on;