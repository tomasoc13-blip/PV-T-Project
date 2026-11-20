clc; close all;

D_i = 18*(10^-3); 
vwind = 22.4*(1000)/(60^2); 
ro = 998; 
V_litre_min = 3; 
mdot = V_litre_min*ro/60000; 
mu = 0.890*10^-3; 
k_water = 0.6; 

params.G = 437; 
params.Tamb = 15 + 273.15; 
params.Tsky = 7 + 273.15; 
params.h_wind = 2.8 + 3*vwind;
params.Tin = 10 + 273.15;
params.mdot = mdot;
params.Cp = 4182;
params.Ac = 2; 
params.alpha = 0.9; 
params.eps_p = 0.9;
params.P = 3.14*(D_i); 
params.L = 5; 
params.Lgap = 0.05;

params.h_gap = 2; 

Re_D = 4*mdot/(3.14*D_i*mu);
Pr = 4182*mu/k_water;



if Re_D < 2300
    Nu_D = 3.66; 
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

Tank_volume = 200/1000;        
M_tank = ro * Tank_volume;     
T_tank = 10 + 273.15;          
target_temp = 45 + 273.15;     

dt = 60;                       
max_hours = 72;                
max_steps = max_hours * 60;    

time_history = zeros(max_steps, 1);
tank_temp_history = zeros(max_steps, 1);
tout_history = zeros(max_steps, 1);
tp_history = zeros(max_steps, 1);     
tg_history = zeros(max_steps, 1);     
h_gap_history = zeros(max_steps, 1);
quseful_history = zeros(max_steps, 1);
qelec_history = zeros(max_steps, 1);
qloss_conv_history = zeros(max_steps, 1);
qloss_rad_history = zeros(max_steps, 1);
eta_th_history = zeros(max_steps, 1);
eta_elec_history = zeros(max_steps, 1);
eta_total_history = zeros(max_steps, 1);

options = optimset('Display', 'off'); 
fprintf('\n--- RUNNING 200L TANK SIMULATION (1-MIN INCREMENTS) ---\n');

for step = 1:max_steps
    params.Tin = T_tank;
    
    Tp_guess = T_tank + 10; 
    Tp_solved = fzero(@(Tp) glazed_residual(Tp, params), Tp_guess, options);
    
    [~, Tout_solved, Q_useful_solved, h_gap_solved, Q_elec_solved, eta_elec, Q_loss_rad, Q_loss_conv, Tg_solved] = glazed_residual(Tp_solved, params);
    

    elapsed_hours = (step - 1) / 60;
    time_history(step) = elapsed_hours;
    tank_temp_history(step) = T_tank - 273.15;
    tout_history(step) = Tout_solved - 273.15;
    tp_history(step) = Tp_solved - 273.15;  
    tg_history(step) = Tg_solved - 273.15;  
    h_gap_history(step) = h_gap_solved; 
    quseful_history(step) = Q_useful_solved;
    qelec_history(step) = Q_elec_solved;
    qloss_conv_history(step) = Q_loss_conv;
    qloss_rad_history(step) = Q_loss_rad;
    
    eta_th = Q_useful_solved / (params.G * params.Ac);
    eta_th_history(step) = eta_th * 100;
    eta_elec_history(step) = eta_elec * 100;
    eta_total_history(step) = eta_th_history(step) + eta_elec_history(step);
    
    if mod(step - 1, 15) == 0
        fprintf('Time: %5.2f hrs | Tank: %5.2f °C | Thermal: %5.1f W | Elec: %5.1f W | Tp: %5.2f °C | Tg: %5.2f °C\n', ...
            elapsed_hours, tank_temp_history(step), Q_useful_solved, Q_elec_solved, tp_history(step), tg_history(step));
    end
    
    if T_tank >= target_temp
        fprintf('\nTarget 45 °C successfully reached at %.2f hours (%d minutes)!\n', elapsed_hours, step - 1);
        time_history = time_history(1:step);
        tank_temp_history = tank_temp_history(1:step);
        tout_history = tout_history(1:step);
        tp_history = tp_history(1:step);
        tg_history = tg_history(1:step);
        quseful_history = quseful_history(1:step);
        qelec_history = qelec_history(1:step);
        qloss_conv_history = qloss_conv_history(1:step);
        qloss_rad_history = qloss_rad_history(1:step);
        eta_th_history = eta_th_history(1:step);
        eta_elec_history = eta_elec_history(1:step);
        eta_total_history = eta_total_history(1:step);
        break;
    end
    
    P_pump = 8; 
    
    Q_total_watts = Q_useful_solved + (Q_elec_solved - P_pump); 
    Q_total_joules = Q_total_watts * dt;
    T_tank = T_tank + Q_total_joules / (M_tank * params.Cp);
end

figure('Color', 'white', 'Position', [100, 100, 950, 600]);
plot(time_history, tank_temp_history, 'r-', 'LineWidth', 2);
hold on;
plot(time_history, tout_history, 'b--', 'LineWidth', 1.5);
plot(time_history, tp_history, 'm-', 'LineWidth', 1.5);
plot(time_history, tg_history, 'g-.', 'LineWidth', 1.5);
yline(45, 'k:', 'LineWidth', 1.5, 'Label', '45°C Target');

xlabel('Time (Hours)', 'FontSize', 12);
ylabel('Temperature (°C)', 'FontSize', 12);
title('PVT System & 200L Tank Transient Temperatures', 'FontSize', 14);
legend('Tank Temp (T_{tank})', 'Collector Outlet (T_{out})', 'Absorber Plate (T_p)', 'Glass Cover (T_g)', 'Target (45°C)', 'Location', 'east');
grid on;

zoom_idx = find(time_history >= 14 & time_history <= 22);
if isempty(zoom_idx)
    zoom_idx = 1:length(time_history);
end

inset_ax = axes('Position', [0.19, 0.55, 0.2, 0.2]);
box(inset_ax, 'on');

plot(inset_ax, time_history(zoom_idx), tank_temp_history(zoom_idx), 'r-', 'LineWidth', 1.5);
hold(inset_ax, 'on');
plot(inset_ax, time_history(zoom_idx), tout_history(zoom_idx), 'b--', 'LineWidth', 1.2);
plot(inset_ax, time_history(zoom_idx), tp_history(zoom_idx), 'm-', 'LineWidth', 1.0);
plot(inset_ax, time_history(zoom_idx), tg_history(zoom_idx), 'g-.', 'LineWidth', 1.0);
yline(inset_ax, 45, 'k:', 'LineWidth', 1.0);

title(inset_ax, 'Zoom: 45°C Target', 'FontSize', 8, 'FontWeight', 'bold');
xlabel(inset_ax, 'Time (hrs)', 'FontSize', 7);
ylabel(inset_ax, 'Temp (°C)', 'FontSize', 7);
grid(inset_ax, 'on');

xlim(inset_ax, [14, 22]);
ylim(inset_ax, [35, 55]);

figure('Color', 'white', 'Position', [150, 150, 950, 700]);

subplot(2, 1, 1);
plot(time_history, quseful_history, 'g-', 'LineWidth', 1.8);
hold on;
plot(time_history, qelec_history, 'c-', 'LineWidth', 1.5);
plot(time_history, qloss_conv_history, 'r--', 'LineWidth', 1.5);
plot(time_history, qloss_rad_history, 'k:', 'LineWidth', 1.5);
xlabel('Time (Hours)', 'FontSize', 11);
ylabel('Power / Heat Loss (W)', 'FontSize', 11);
title('Transient Power Generation and Thermal Losses', 'FontSize', 13, 'FontWeight', 'bold');
legend('Useful Thermal (Q_{useful})', 'Electrical (Q_{elec})', 'Convective Loss (Q_{loss,conv})', 'Radiative Loss (Q_{loss,rad})', 'Location', 'north', 'FontSize', 9);
grid on;
box off;

subplot(2, 1, 2);
plot(time_history, eta_th_history, 'r-', 'LineWidth', 1.8);
hold on;
plot(time_history, eta_elec_history, 'b-', 'LineWidth', 1.5);
plot(time_history, eta_total_history, 'k-', 'LineWidth', 2);
mean_eta_total = mean(eta_total_history);
yline(mean_eta_total, 'k--', 'LineWidth', 1, 'DisplayName', sprintf('Avg Combined Efficiency (%.1f%%)', mean_eta_total));
xlabel('Time (Hours)', 'FontSize', 11);
ylabel('Efficiency (%)', 'FontSize', 11);
title('PVT Collector Instantaneous Efficiencies Over Time', 'FontSize', 13, 'FontWeight', 'bold');
legend('Thermal Efficiency (\eta_{th})', 'Electrical Efficiency (\eta_{elec})', 'Combined Efficiency (\eta_{total})', sprintf('Time Avg Combined Efficiency (%.1f%%)', mean_eta_total), 'Location', 'northeast', 'FontSize', 9);
grid on;
box off;

% Lgap_vec = linspace(0.004, 0.1, 40);      
% avg_Q_loss_conv = NaN(size(Lgap_vec));
% heating_time_Lgap = NaN(size(Lgap_vec));
% avg_eta_total = NaN(size(Lgap_vec));
% P_pump = 8;
% 
% fprintf('\n--- RUNNING TRANSIENT L_gap SWEEP (4 mm to 100 mm) ---\n');
% 
% for k = 1:length(Lgap_vec)
%     params_sweep = params;
%     params_sweep.Lgap = Lgap_vec(k);
% 
%     T_tank_sweep = 10 + 273.15;
%     sum_Q_conv = 0;
%     sum_eta_total = 0;
% 
%     for step = 1:max_steps
%         params_sweep.Tin = T_tank_sweep;
%         Tp_guess = T_tank_sweep + 10;
% 
%         Tp_solved = fzero(@(Tp) glazed_residual(Tp, params_sweep), Tp_guess, options);
%         [~, ~, Q_useful_solved, ~, Q_elec_solved, eta_elec, ~, Q_loss_conv_solved] = glazed_residual(Tp_solved, params_sweep);
% 
%         sum_Q_conv = sum_Q_conv + Q_loss_conv_solved;
%         eta_th = Q_useful_solved / (params_sweep.G * params_sweep.Ac);
%         sum_eta_total = sum_eta_total + (eta_th + eta_elec) * 100;
% 
%         if T_tank_sweep >= target_temp
%             avg_Q_loss_conv(k)   = sum_Q_conv / step;
%             heating_time_Lgap(k) = (step - 1) / 60;
%             avg_eta_total(k)     = sum_eta_total / step;
%             break;
%         end
% 
%         Q_total_watts = Q_useful_solved + (Q_elec_solved - P_pump);
%         T_tank_sweep = T_tank_sweep + (Q_total_watts * dt) / (M_tank * params_sweep.Cp);
%     end
% end
% 
% figure('Color', 'white', 'Position', [100, 100, 1300, 450]);
% 
% subplot(1, 2, 1);
% plot(Lgap_vec * 1000, avg_Q_loss_conv, 'r-', 'LineWidth', 2);
% hold on;
% xline(50, 'k--', 'LineWidth', 1.5, 'Label', '50 mm Design Point');
% xlabel('Air Gap Thickness, L_{gap} (mm)', 'FontSize', 11, 'FontWeight', 'bold');
% ylabel('Time-Averaged Q_{loss,conv} (W)', 'FontSize', 11, 'FontWeight', 'bold');
% title('Average Convective Heat Loss', 'FontSize', 12);
% grid on;
% 
% subplot(1, 2, 2);
% plot(Lgap_vec * 1000, heating_time_Lgap, 'b-', 'LineWidth', 2);
% hold on;
% xline(50, 'k--', 'LineWidth', 1.5, 'Label', '50 mm Design Point');
% xlabel('Air Gap Thickness, L_{gap} (mm)', 'FontSize', 11, 'FontWeight', 'bold');
% ylabel('Time to Reach 45°C (Hours)', 'FontSize', 11, 'FontWeight', 'bold');
% title('Tank Warm-Up Time', 'FontSize', 12);
% grid on;



% Plot Losses vs Tp-


qloss_total_history = qloss_conv_history + qloss_rad_history;

figure('Color', 'white', 'Position', [180, 120, 980, 780]);

%% --- SUBPLOT 1: HEAT LOSSES VS. Tp (TOTAL COLOUR-CODED BY TIME) ---
subplot(2, 1, 1);
hold on;

plot(tp_history, qloss_conv_history, 'r--', 'LineWidth', 1.8, 'DisplayName', 'Plate Convective Loss (Q_{loss,conv})');
plot(tp_history, qloss_rad_history,  'b:',  'LineWidth', 1.8, 'DisplayName', 'Plate Radiative Loss (Q_{loss,rad})');
plot(tp_history, qloss_total_history, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0, 'HandleVisibility', 'off');

scatter(tp_history, qloss_total_history, 30, time_history, 'filled', 'DisplayName', 'Total Loss (colour-coded by time)');

mark_idx = find(mod(round(time_history * 100), 300) == 0);
if mark_idx(end) ~= length(time_history)
    mark_idx(end+1) = length(time_history);
end
plot(tp_history(mark_idx), qloss_total_history(mark_idx), 'ko', 'MarkerSize', 6, 'LineWidth', 1.2, 'HandleVisibility', 'off');

y_span = max(qloss_total_history) - min(qloss_total_history);
y_offset = 0.025 * max(y_span, 1);

for k = 1:length(mark_idx)
    idx = mark_idx(k);
    if k == length(mark_idx)
        lbl = sprintf('%.2fh', time_history(idx));
    else
        lbl = sprintf('%.0fh', time_history(idx));
    end
    text(tp_history(idx), qloss_total_history(idx) + y_offset, lbl, ...
         'FontSize', 9, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

caxis([min(time_history), max(time_history)]); % Use clim(...) in R2022a+ if preferred
xlabel('Absorber Plate Temperature, T_p (°C)', 'FontSize', 11);
ylabel('Heat Loss (W)', 'FontSize', 11);
title('Plate Heat Losses vs. T_p (Total Loss Colour-Coded by Time)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 9);
grid on;
box off;

%% --- SUBPLOT 2: EFFICIENCIES VS. Tp (COMBINED COLOUR-CODED BY TIME) ---
subplot(2, 1, 2);
hold on;

plot(tp_history, eta_th_history,   'r--', 'LineWidth', 1.8, 'DisplayName', 'Thermal Efficiency (\eta_{th})');
plot(tp_history, eta_elec_history, 'b:',  'LineWidth', 1.8, 'DisplayName', 'Electrical Efficiency (\eta_{elec})');
plot(tp_history, eta_total_history, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0, 'HandleVisibility', 'off');

scatter(tp_history, eta_total_history, 30, time_history, 'filled', 'DisplayName', 'Combined Efficiency (colour-coded by time)');

plot(tp_history(mark_idx), eta_total_history(mark_idx), 'ko', 'MarkerSize', 6, 'LineWidth', 1.2, 'HandleVisibility', 'off');

eta_span = max(eta_total_history) - min(eta_total_history);
eta_offset = 0.025 * max(eta_span, 1);

for k = 1:length(mark_idx)
    idx = mark_idx(k);
    if k == length(mark_idx)
        lbl = sprintf('%.2fh', time_history(idx));
    else
        lbl = sprintf('%.0fh', time_history(idx));
    end
    text(tp_history(idx), eta_total_history(idx) + eta_offset, lbl, ...
         'FontSize', 9, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

caxis([min(time_history), max(time_history)]);
xlabel('Absorber Plate Temperature, T_p (°C)', 'FontSize', 11);
ylabel('Efficiency (%)', 'FontSize', 11);
title('PVT Efficiencies vs. T_p (Combined Efficiency Colour-Coded by Time)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 9);
grid on;
box off;

%% --- SHARED THIN COLORBAR ON THE RIGHT ---
cb = colorbar('Position', [0.92, 0.14, 0.012, 0.74]);
cb.Label.String = 'Elapsed Time (Hours)';
cb.Label.FontSize = 11;