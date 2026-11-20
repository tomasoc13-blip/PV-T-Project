clc; clear; close all;

%  Define Parametric Sweep Ranges ---
flow_vec   = 0.1:0.1:3.9;                      % 25 elements: 0.1 to 4.9 in 0.2 steps (5.0 isn't hit cleanly from 0.1)
length_vec = linspace(1, 12, length(flow_vec)); % 25 elements: evenly spaced from 1 to 10 m

nF = length(flow_vec);
nL = length(length_vec);


Time_Matrix = zeros(nL, nF);

fprintf('Running %d simulations (Flow Rate x Pipe Length)...\n', nF * nL);

%  Nested Parametric Loop ---
for i = 1:nL
    for j = 1:nF
        L_val = length_vec(i);
        V_val = flow_vec(j);
        
        Time_Matrix(i, j) = sim_tank_parametric(V_val, L_val, 20 * (10^-3));
    end
end


[FlowGrid, LengthGrid] = meshgrid(flow_vec, length_vec);

figure('Color', 'white', 'Position', [100, 100, 1200, 500]);

% 3D Parametric Surface
subplot(1, 1, 1);
surf(FlowGrid, LengthGrid, Time_Matrix, 'EdgeColor', 'k', 'FaceAlpha', 0.9);
colormap(jet);
cb1 = colorbar;
cb1.Label.String = 'Warm-Up Time (Hours)';
xlabel('Flow Rate, V (L/min)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Pipe Length, L (m)', 'FontSize', 14, 'FontWeight', 'bold');
zlabel('Time to 45°C (Hours)', 'FontSize', 14, 'FontWeight', 'bold');
title('3D Parametric Surface (L_{gap} = 0.03 m, P_{pump} = 8 W)', 'FontSize', 16);
view(-35, 30);
grid on;



fprintf('\n--- WARM-UP TIME MATRIX (HOURS) ---\n');
ResultTable = array2table(Time_Matrix, ...
    'RowNames', compose('L = %d m', length_vec), ...
    'VariableNames', compose('V_%.1f_Lmin', flow_vec));
disp(ResultTable);

% 
% SAME SCRIPT BELOW AS BEFORE, just copied it here for parametric sweep
% 
% 
function hours_to_45 = sim_tank_parametric(V_litre_min, L_pipe, D_i)
    %D_i = 20 * (10^-3);
    vwind = 31 * (1000) / (60^2);
    ro = 998;
    mu = 0.890 * 10^-3;
    k_water = 0.6;
    mdot = V_litre_min * ro / 60000;

    params.G = 437;
    params.Tamb = 15 + 273.15;
    params.Tsky = 7 + 273.15;
    params.h_wind = 2.8 + 3 * vwind;
    params.Tin = 10 + 273.15;
    params.mdot = mdot;
    params.Cp = 4182;
    params.Ac = 2;
    params.alpha = 0.9;
    params.eps_p = 0.9;
    params.P = 3.14 * (D_i);
    params.L = L_pipe;
    params.Lgap = 0.02;
    

    % Internal tube flow calculations (Gnielinski / Laminar)
    Re_D = 4 * mdot / (3.14 * D_i * mu);
    Pr = 4182 * mu / k_water;
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
    params.hfi = k_water * Nu_D / (D_i);

    % Closed-loop 200L tank simulation (1-min step)
    Tank_volume = 200 / 1000;
    M_tank = ro * Tank_volume;
    T_tank = 10 + 273.15;
    target_temp = 45 + 273.15;
    dt = 60;
    max_hours = 72;
    max_steps = max_hours * 60;
    options = optimset('Display', 'off');
    P_pump = 8;

    hours_to_45 = NaN;

    for step = 1:max_steps
        params.Tin = T_tank;
        Tp_guess = T_tank + 10;
        Tp_solved = fzero(@(Tp) glazed_residual(Tp, params), Tp_guess, options);
        [~, ~, Q_useful_solved, ~, Q_elec_solved] = glazed_residual(Tp_solved, params);

        if T_tank >= target_temp
            hours_to_45 = (step - 1) / 60;
            return;
        end

        Q_total_watts = Q_useful_solved + (Q_elec_solved - P_pump);
        Q_total_joules = Q_total_watts * dt;
        T_tank = T_tank + Q_total_joules / (M_tank * params.Cp);
    end
end





% Diameter Sweep

V_fixed = 3.0; % L/min
L_fixed = 5.0; % m

% 20 points from 5 mm to 23 mm, and 50 points from 23 mm to 40 mm to
% capture the transitional points correctly
d_low  = linspace(0.005, 0.023, 20);
d_high = linspace(0.023, 0.040, 50);

diam_sweep_vec = [d_low, d_high(2:end)];

time_diam_sweep = zeros(size(diam_sweep_vec));

fprintf('\n--- RUNNING DIAMETER SENSITIVITY SWEEP (V = %.1f L/min, L = %.1f m) ---\n', V_fixed, L_fixed);

for d = 1:length(diam_sweep_vec)
    D_val = diam_sweep_vec(d);
    time_diam_sweep(d) = sim_tank_parametric(V_fixed, L_fixed, D_val);
end




figure('Color', 'white', 'Position', [150, 150, 850, 480]);
plot(diam_sweep_vec * 1000, time_diam_sweep, 'b-', 'LineWidth', 2);
hold on;

xline(23.8, '--r', 'LineWidth', 1.5);
xline(31.16, '--k', 'LineWidth', 1.5);


yl = ylim;
y_top = yl(1) + 0.90 * (yl(2) - yl(1));
y_mid = yl(1) + 0.50 * (yl(2) - yl(1));


text(23.4, y_top-0.5, 'Re = 3000 (Turbulent)', 'Color', 'r', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'right', 'FontSize', 11);
text(31.4+1, y_top-0.5, 'Re = 2300 (Laminar)', 'Color', 'k', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'FontSize', 11);


text(27.4, y_mid+0.5, 'Transitional (Interpolation)', 'Color', [0.3 0.3 0.3], 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90, 'FontSize', 12);

xlabel('Inner Tube Diameter, D_i (mm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Time to Reach 45°C (Hours)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Warm-Up Time vs. Tube Diameter (at V = %.1f L/min, L = %.1f m)', V_fixed, L_fixed), 'FontSize', 13);
grid off;
hold off;