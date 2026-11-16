clc, close all, clear

% Define inputs
D_i = 20 * (10^-3); % Inner diameter (m)
vwind = 20 * (1000) / (60^2); % Wind speed (m/s)
ro = 998; % Water density (kg/m^3)
V_litre_min = 2; % Flow rate (L/min)
mdot = V_litre_min * ro / 60000; % Mass flow rate (kg/s)
mu = 0.890 * 10^-3; % Dynamic viscosity (Pa.s)
k_water = 0.6; % Thermal conductivity (W/m.K)

% Collector & ambient parameters
params.G = 436; 
params.Tamb = 17 + 273.15; 
params.Tsky = 15 + 273.15; 
params.h_wind = 2.8 + 3 * vwind;
params.mdot = mdot;
params.Cp = 4182;
params.Ac = 2;
params.alpha = 0.9; % Kirchhoff's Law
params.eps_p = 0.9;
params.P = 3.14 * (D_i); % Wetted perimeter (m)
params.L = 10; % Length (m)

% Internal convection coefficient (hfi)
Re_D = 4 * mdot / (3.14 * D_i * mu);
Pr = 4182 * mu / k_water;
if Re_D < 2300
    Nu_D = 3.66; % Fully developed laminar
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
params.hfi = k_water * Nu_D / (D_i);

% Tank & Euler Time-Stepping Setup 
V_tank_L = 200; % Tank volume (Litres)
m_tank = (V_tank_L / 1000) * ro; % Tank water mass (kg)
T_tank = 10 + 273.15; % Initial tank temperature (K)
T_target = 45 + 273.15; % Target temperature (45 °C in K)
dt = 3600; % Time step: 1 hour (seconds)
time_hrs = 0;
max_hrs = 200; 

options = optimset('Display', 'off'); 
Tp_guess = 310;

% Arrays
hist_time = [0];
hist_Ttank = [T_tank - 273.15];
hist_Tplate = [];
hist_Quseful = [];
hist_Qloss_conv = [];
hist_Qloss_rad = [];
hist_Qelec = [];

fprintf('Time (hr) | T_tank (°C) | T_plate (°C) | Q_useful (W) | Q_elec (W) | Total Power (W) | Q_loss_conv W | Q_loss_rad W\n');
fprintf('-------------------------------------------------------------------------------------\n');

while T_tank < T_target && time_hrs < max_hrs
    params.Tin = T_tank;
    
    Tp_solved = fzero(@(Tp) bare_residual(Tp, params), Tp_guess);
    [~, Tout_solved, Q_useful_solved, Q_elec_solved, Q_loss_conv, Q_loss_rad] = bare_residual(Tp_solved, params);
    
    Tp_guess = Tp_solved;
    
    Q_thermal_net = max(0, Q_useful_solved);
    if Q_thermal_net == 0 
        params.hfi = 0;
    end
    Q_total = Q_thermal_net + Q_elec_solved;
    
    fprintf('%9.1f | %11.2f | %12.2f | %12.2f | %10.2f | %15.2f | %10.2f | %10.2f\n', ...
        time_hrs, T_tank - 273.15, Tp_solved - 273.15, Q_thermal_net, Q_elec_solved, Q_total, Q_loss_conv, Q_loss_rad);
    
    % Store history
    hist_Tplate(end+1) = Tp_solved - 273.15;
    hist_Quseful(end+1) = Q_thermal_net;
    hist_Qloss_conv(end+1) = Q_loss_conv;
    hist_Qloss_rad(end+1) = Q_loss_rad;
    hist_Qelec(end+1) = Q_elec_solved;
    
    dT = (Q_total * dt) / (m_tank * params.Cp);
    T_tank = T_tank + dT;
    time_hrs = time_hrs + (dt / 3600);
    
    hist_time(end+1) = time_hrs;
    hist_Ttank(end+1) = T_tank - 273.15;
end
% Pad Tplate array to match length of hist_time for plotting
hist_Tplate(end+1) = hist_Tplate(end);

fprintf('-------------------------------------------------------------------------------------\n');
fprintf('Final Tank Temp: %.2f °C reached at %.1f hours\n', T_tank - 273.15, time_hrs);

%
figure('Position', [100, 100, 800, 700]);

% Subplot 1: Temperatures
subplot(2,1,1);
plot(hist_time, hist_Ttank, 'b-', 'LineWidth', 2, 'DisplayName', 'Tank Temp (T_{tank})');
hold on;
plot(hist_time, hist_Tplate, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Plate Temp (T_{plate})');
yline(45, 'r:', 'Target 45 °C', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(9, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(9.2, 25, 'Pump shutoff', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
ylabel('Temperature (°C)');
title('PVT System Temperatures');
legend('Location', 'southeast', 'FontSize', 10);
grid on;
box off;

% Subplot 2: Power and Losses
subplot(2,1,2);
plot(hist_time(1:end-1), hist_Quseful, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Useful Thermal Gain (Q_{useful})');
hold on;
plot(hist_time(1:end-1), hist_Qelec, 'c-', 'LineWidth', 1.5, 'DisplayName', 'Electrical Output (Q_{elec})');
plot(hist_time(1:end-1), hist_Qloss_conv, 'm-.', 'LineWidth', 1.5, 'DisplayName', 'Convective Loss (Q_{loss,conv})');
plot(hist_time(1:end-1), hist_Qloss_rad, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Radiative Loss (Q_{loss,rad})');
xline(9, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(9.2, 250, 'Pump shutoff', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
xlabel('Time (hours)');
ylabel('Power / Heat Transfer (W)');
title('Power Generation and Thermal Losses');
legend('Location', 'northeast', 'FontSize', 10);
grid on;
box off;