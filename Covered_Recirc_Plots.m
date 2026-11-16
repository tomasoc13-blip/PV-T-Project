clc; close all;

% Define inputs
D_i = 20*(10^-3); % Convert mm to m, inner diameter
vwind = 20*(1000)/(60^2); % convert kmh to m/s
ro = 998; % water density, kg/m^3
V_litre_min = 2; % Design parameter
mdot = V_litre_min*ro/60000; % kg/sec
mu = 0.890*10^-3; % Dynamic viscosity of water
k_water = 0.6; % Thermal conductivity of water @ approx 25 celsius

params.G = 436; % Solar irradiation 
params.Tamb = 15 + 273.15; 
params.Tsky = 10 + 273.15; 
params.h_wind = 2.8 + 3*vwind;
params.Tin = 10 + 273.15;
params.mdot = mdot;
params.Cp = 4182;
params.Ac = 2; % Assigned absorber area
params.alpha = 0.9; 
params.eps_p = 0.9;
params.P = 3.14*(D_i); % Wetted perimeter
params.L = 10; % Total length of piping
params.Lgap = 0.02; % Air gap thickness between plate and glass


% Internal tube flow calculations (Gnielinski / Laminar)
Re_D = 4*mdot/(3.14*D_i*mu);
Pr = 4182*mu/k_water;
if Re_D < 2300
    Nu_D = 3.66; % Const wall temp
elseif Re_D < 3000
    f_3000 = (0.790 * log(3000) - 1.64)^(-2);
    Nu_3000 = ((f_3000 / 8) * (3000 - 1000) * Pr) / ...
              (1 + 12.7 * sqrt(f_3000 / 8) * (Pr^(2/3) - 1));
    gamma = (Re_D - 2300) / (3000 - 2300);
    Nu_D = (1 - gamma) * 4.36 + gamma * Nu_3000;
else
    f = (0.790 * log(Re_D) - 1.64)^(-2);
    Nu_D = ((f / 8) * (Re_D - 1000) * Pr) / ...
         (1 + 12.7 * sqrt(f / 8) * (Pr^(2/3) - 1));
end
params.hfi = k_water*Nu_D/(D_i);

% Tank Recirculation Simulation
Tank_volume = 200/1000;        % 200 litres = 0.2 m^3
M_tank = ro * Tank_volume;     % Mass of water in the tank (kg)
T_tank = 10 + 273.15;          % Initial tank temperature (K)
target_temp = 45 + 273.15;     % Target temperature (K)
dt = 60;                       % Time step in seconds (1 minute)
max_hours = 72;                % Max simulation window (hours)
max_steps = max_hours * 60;    % Total 1-minute steps

%  history arrays
time_history = zeros(max_steps, 1);
tank_temp_history = zeros(max_steps, 1);
tout_history = zeros(max_steps, 1);
tplate_history = zeros(max_steps, 1);
quseful_history = zeros(max_steps, 1);
qelec_history = zeros(max_steps, 1);
qloss_conv_history = zeros(max_steps, 1);
qloss_rad_history = zeros(max_steps, 1);
h_gap_history = zeros(max_steps, 1);

options = optimset('Display', 'off'); 
fprintf('\n--- RUNNING 200L TANK SIMULATION (1-MIN INCREMENTS) ---\n');

for step = 1:max_steps
    % 1. Inlet water comes from the current tank temperature
    params.Tin = T_tank;
    
    % 2. Solve for plate temperature at this minute's conditions
    Tp_guess = T_tank + 10; 
    Tp_solved = fzero(@(Tp) glazed_residual(Tp, params), Tp_guess, options);
    
    % 3. Extract outputs for this minute
    [~, Tout_solved, Q_useful_solved, h_gap_solved, Q_elec_solved, eta_elec, Q_loss_rad, Q_loss_conv] = glazed_residual(Tp_solved, params);
    
    % Store history
    elapsed_hours = (step - 1) / 60;
    time_history(step) = elapsed_hours;
    tank_temp_history(step) = T_tank - 273.15;
    tout_history(step) = Tout_solved - 273.15;
    tplate_history(step) = Tp_solved - 273.15;
    quseful_history(step) = Q_useful_solved;
    qelec_history(step) = Q_elec_solved;
    qloss_conv_history(step) = Q_loss_conv;
    qloss_rad_history(step) = Q_loss_rad;
    h_gap_history(step) = h_gap_solved;
    
    % Print progress
    if mod(step - 1, 15) == 0
        fprintf('Time: %5.2f hrs | Tank: %5.2f °C | Thermal: %5.1f W | Elec: %5.1f W | T_p: %5.2f °C\n', ...
            elapsed_hours, tank_temp_history(step), Q_useful_solved, Q_elec_solved, tplate_history(step));
    end
    
    % Check if we hit our 45 °C target
    if T_tank >= target_temp
        fprintf('\nTarget 45 °C successfully reached at %.2f hours (%d minutes)!\n', elapsed_hours, step - 1);
        % Trim history vectors to actual completion step
        time_history = time_history(1:step);
        tank_temp_history = tank_temp_history(1:step);
        tout_history = tout_history(1:step);
        tplate_history = tplate_history(1:step);
        quseful_history = quseful_history(1:step);
        qelec_history = qelec_history(1:step);
        qloss_conv_history = qloss_conv_history(1:step);
        qloss_rad_history = qloss_rad_history(1:step);
        h_gap_history = h_gap_history(1:step);
        break;
    end
    
    % 4. Update tank temperature using thermal heat + electrical immersion boost
    P_pump = 8; % Pump power, assumed constant draw
    Q_total_watts = Q_useful_solved + (Q_elec_solved - P_pump); 
    Q_total_joules = Q_total_watts * dt;
    T_tank = T_tank + Q_total_joules / (M_tank * params.Cp);
end

% --- PLOT THE RESULTS WITH SUBPLOTS ---
figure('Color', 'white', 'Position', [100, 100, 900, 700]);

% Subplot 1: Temperatures
subplot(2,1,1);
plot(time_history, tank_temp_history, 'r-', 'LineWidth', 2, 'DisplayName', 'Tank Temp (T_{tank})');
hold on;

plot(time_history, tout_history, 'b:', 'LineWidth', 1.5, 'DisplayName', 'Outlet Temp (T_{out})');
yline(45, 'k:', 'LineWidth', 1.5, 'Label', '45°C Target', 'HandleVisibility', 'off');
ylabel('Temperature (°C)', 'FontSize', 11);
title('PVT System Temperatures', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 10);
grid on;
box off;

% Subplot 2: Power and Losses
subplot(2,1,2);
plot(time_history, quseful_history, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Useful Thermal (Q_{useful})');
hold on;
plot(time_history, qelec_history, 'c-', 'LineWidth', 1.5, 'DisplayName', 'Electrical Output (Q_{elec})');
plot(time_history, qloss_conv_history, 'r-.', 'LineWidth', 1.5, 'DisplayName', 'Convective Loss (Q_{loss,conv})');
plot(time_history, qloss_rad_history, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Radiative Loss (Q_{loss,rad})');
xlabel('Time (Hours)', 'FontSize', 11);
ylabel('Power / Heat Loss (W)', 'FontSize', 11);
title('Collector Power and Thermal Losses', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'north', 'FontSize', 10);
grid on;
box off;