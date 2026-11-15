clc; close all;

% Define inputs
D_i = 20*(10^-3); % Convert mm to m, inner diameter
vwind = 20*(1000)/(60^2); % convert kmh to m/s
ro = 998; % water density, kg/m^3
V_litre_min = 3; % Design parameter
mdot = V_litre_min*ro/60000; % kg/sec
mu = 0.890*10^-3; % Dynamic viscosity of water
k_water = 0.6; % Thermal conductivity of water @ approx 25 celsius

params.G = 436;%437; % Solar irradiation 
params.Tamb = 17 + 273.15; 
params.Tsky = 15 + 273.15; 
params.h_wind = 2.8 + 3*vwind;
params.Tin = 10 + 273.15;
params.mdot = mdot;
params.Cp = 4182;
params.Ac = 2; % Assigned absorber area
params.alpha = 0.9; 
params.eps_p = 0.9;
params.P = 3.14*(D_i); % Wetted perimeter
params.L = 10; % Total length of piping

% 
params.h_gap = 2; % Convective heat transfer coefficient for the stagnant air gap (W/m^2K)

% Internal tube flow calculations (Gnielinski / Laminar)
Re_D = 4*mdot/(3.14*D_i*mu);
Pr = 4182*mu/k_water;

if Re_D < 2300
    Nu_D = 3.66; % Const wall temp

elseif Re_D < 3000
    % Gnielinski evaluated at Re = 3000
    f_3000 = (0.790 * log(3000) - 1.64)^(-2);
    Nu_3000 = ((f_3000 / 8) * (3000 - 1000) * Pr) / ...
              (1 + 12.7 * sqrt(f_3000 / 8) * (Pr^(2/3) - 1));
    
    % Linear interpolation across 2300 <= Re < 3000
    gamma = (Re_D - 2300) / (3000 - 2300);
    Nu_D = (1 - gamma) * 4.36 + gamma * Nu_3000;
else
    % Standard Gnielinski for Re >= 3000
    f = (0.790 * log(Re_D) - 1.64)^(-2);
    Nu_D = ((f / 8) * (Re_D - 1000) * Pr) / ...
         (1 + 12.7 * sqrt(f / 8) * (Pr^(2/3) - 1));
end

params.hfi = k_water*Nu_D/(D_i);


% --- CLOSED-LOOP TANK SIMULATION (1-MINUTE TIME STEP) ---
Tank_volume = 200/1000;        % 200 litres = 0.2 m^3, fixed design variable
M_tank = ro * Tank_volume;     % Mass of water in the tank (kg)
T_tank = 10 + 273.15;          % Initial tank temperature (K)
target_temp = 45 + 273.15;     % Target temperature (K)

dt = 60;                       % Time step in seconds (1 minute)
max_hours = 72;                % Max simulation window (hours)
max_steps = max_hours * 60;    % Total 1-minute steps (1440)

time_history = zeros(max_steps, 1);
tank_temp_history = zeros(max_steps, 1);
tout_history = zeros(max_steps, 1);

options = optimset('Display', 'off'); 
fprintf('\n--- RUNNING 200L TANK SIMULATION (1-MIN INCREMENTS) ---\n');

for step = 1:max_steps
    % 1. Inlet water comes from the current tank temperature
    params.Tin = T_tank;
    
    % 2. Solve for plate temperature at this minute's conditions
    Tp_guess = T_tank + 10; 
    Tp_solved = fzero(@(Tp) glazed_residual(Tp, params), Tp_guess, options);
    
% 3. Extract outputs for this minute (now including Q_elec)
    [~, Tout_solved, Q_useful_solved, h_gap_solved, Q_elec_solved, eta_elec] = glazed_residual(Tp_solved, params);
    
    % Store history 
    elapsed_hours = (step - 1) / 60;
    time_history(step) = elapsed_hours;
    tank_temp_history(step) = T_tank - 273.15;
    tout_history(step) = Tout_solved - 273.15;
    
    % Print progress including electrical power
    if mod(step - 1, 15) == 0
        fprintf('Time: %5.2f hrs | Tank: %5.2f °C | Thermal: %5.1f W | Elec: %5.1f W | Eta_elec: %5.3f\n', ...
            elapsed_hours, tank_temp_history(step), Q_useful_solved, Q_elec_solved, eta_elec);
    end
    
    % Check if we hit our 45 °C target
    if T_tank >= target_temp
        fprintf('\nTarget 45 °C successfully reached at %.2f hours (%d minutes)!\n', elapsed_hours, step - 1);
        time_history = time_history(1:step);
        tank_temp_history = tank_temp_history(1:step);
        tout_history = tout_history(1:step);
        break;
    end
    
    % 4. Update tank temperature using BOTH thermal heat and electric immersion boost
    % (Assume 99% efficiency for converting electricity via an immersion heater)
    % or Assume Heat pump, COP = 3.5
    COP = 3.5;
    Q_total_watts = Q_useful_solved + (COP * Q_elec_solved); 
    Q_total_joules = Q_total_watts * dt;
    T_tank = T_tank + Q_total_joules / (M_tank * params.Cp);
end

% --- PLOT THE RESULTS ---
figure('Color', 'white');
plot(time_history, tank_temp_history, 'r-', 'LineWidth', 2);
hold on;
plot(time_history, tout_history, 'b--', 'LineWidth', 1.5);
yline(45, 'k:', 'LineWidth', 1.5, 'Label', '45°C Target');
xlabel('Time (Hours)', 'FontSize', 12);
ylabel('Temperature (°C)', 'FontSize', 12);
title('200L Tank Heat-Up Curve (1-Minute Resolution)', 'FontSize', 14);
legend('Tank Temp (T_{tank})', 'Collector Outlet (T_{out})', 'Target', 'Location', 'southeast');
grid on;