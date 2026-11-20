clc; close all; clear;

D_i_g = 18 * (10^-3);
L_g   = 5;
vwind_g = 22.4 * (1000) / (60^2);
ro = 998;
V_litre_min_g = 3;
mdot_g = V_litre_min_g * ro / 60000;
mu = 0.890 * 10^-3;
k_water = 0.6;

params_pvt.G = 437;
params_pvt.Tamb = 15 + 273.15;
params_pvt.Tsky = 7 + 273.15;
params_pvt.h_wind = 2.8 + 3 * vwind_g;
params_pvt.Tin = 10 + 273.15;
params_pvt.mdot = mdot_g;
params_pvt.Cp = 4182;
params_pvt.Ac = 2;
params_pvt.alpha = 0.9;
params_pvt.eps_p = 0.9;
params_pvt.P = 3.14 * D_i_g;
params_pvt.L = L_g;
params_pvt.Lgap = 0.05;
params_pvt.run_elec = true;
params_pvt.hfi = calc_hfi(mdot_g, D_i_g, mu, k_water, 4182);

Tank_volume = 200 / 1000;
M_tank = ro * Tank_volume;
target_temp = 45 + 273.15;

dt_g = 60;
max_hours_g = 72;
max_steps_g = max_hours_g * 60;
options = optimset('Display', 'off');
P_pump = 8;

[m_pvt, traj_pvt] = run_glazed_metrics(params_pvt, M_tank, target_temp, dt_g, max_steps_g, P_pump, options);

params_th = params_pvt;
params_th.run_elec = false;
params_th.eps_p = 0.9;
[m_th, traj_th] = run_glazed_metrics(params_th, M_tank, target_temp, dt_g, max_steps_g, P_pump, options);

params_sel = params_pvt;
params_sel.run_elec = false;
params_sel.eps_p = 0.3;
[m_sel, traj_sel] = run_glazed_metrics(params_sel, M_tank, target_temp, dt_g, max_steps_g, P_pump, options);

D_i_u = D_i_g;
L_u   = L_g;
vwind_u = 20 * (1000) / (60^2);
V_litre_min_u = 2;
mdot_u = V_litre_min_u * ro / 60000;

params_ung.G = 436;
params_ung.Tamb = 17 + 273.15;
params_ung.Tsky = 15 + 273.15;
params_ung.h_wind = 2.8 + 3 * vwind_u;
params_ung.mdot = mdot_u;
params_ung.Cp = 4182;
params_ung.Ac = 2;
params_ung.alpha = 0.9;
params_ung.eps_p = 0.9;
params_ung.P = 3.14 * D_i_u;
params_ung.L = L_u;
params_ung.hfi = calc_hfi(mdot_u, D_i_u, mu, k_water, 4182);

T_tank_u = 10 + 273.15;
dt_u = 3600;
time_hrs_u = 0;
max_hrs_u = 200;
t_limit_u = 40;
Tp_guess_u = 310;

eth_u = []; eel_u = []; qloss_u = []; tp_u_vec = []; t_u_vec = [];

while T_tank_u < target_temp && time_hrs_u < max_hrs_u
    if time_hrs_u > t_limit_u
        break;
    end
    params_ung.Tin = T_tank_u;
    Tp_s = fzero(@(Tp) bare_residual(Tp, params_ung), Tp_guess_u, options);
    [~, ~, Q_u_s, Q_e_s, Q_c_s, Q_r_s] = bare_residual(Tp_s, params_ung);
    Tp_guess_u = Tp_s;

    Q_th_net = max(0, Q_u_s);

    t_u_vec(end+1, 1)  = time_hrs_u;
    tp_u_vec(end+1, 1) = Tp_s - 273.15;
    eth_u(end+1, 1)    = (Q_th_net / (params_ung.G * params_ung.Ac)) * 100;
    eel_u(end+1, 1)    = (Q_e_s / (params_ung.G * params_ung.Ac)) * 100;
    qloss_u(end+1, 1)  = Q_c_s + Q_r_s;

    T_tank_u = T_tank_u + ((Q_th_net + Q_e_s) * dt_u) / (M_tank * params_ung.Cp);
    time_hrs_u = time_hrs_u + (dt_u / 3600);
end

m_ung.eta_th = mean(eth_u);
m_ung.eta_el = mean(eel_u);
m_ung.qloss  = mean(qloss_u);
m_ung.time_h = time_hrs_u;

top_lines = {'Unglazed PV', 'Glazed PVT', 'Pure Thermal', 'Pure Thermal'};

eta_th_all  = [m_ung.eta_th, m_pvt.eta_th, m_th.eta_th, m_sel.eta_th];
eta_el_all  = [m_ung.eta_el, m_pvt.eta_el, m_th.eta_el, m_sel.eta_el];
eta_tot_all = eta_th_all + eta_el_all;
time_all    = [m_ung.time_h, m_pvt.time_h, m_th.time_h, m_sel.time_h];
qloss_all   = [m_ung.qloss,  m_pvt.qloss,  m_th.qloss,  m_sel.qloss];

bar_cols = [ ...
    0.18, 0.63, 0.33; ...
    0.12, 0.47, 0.71; ...
    0.84, 0.24, 0.24; ...
    0.58, 0.20, 0.65];

fig = figure('Color', [1, 1, 1], 'InvertHardcopy', 'off', 'Position', [80, 80, 1360, 520]);

ax1 = subplot(1, 3, 1);
set(ax1, 'Color', [1, 1, 1]);
b1 = bar(ax1, [eta_th_all(:), eta_el_all(:)], 'stacked', 'BarWidth', 0.48);
b1(1).FaceColor = [0.85, 0.37, 0.00];
b1(2).FaceColor = [0.12, 0.47, 0.71];
format_case_axis(ax1, top_lines, max(eta_tot_all) * 1.18);
ylabel(ax1, 'Mean Efficiency (%)', 'FontWeight', 'bold');
title(ax1, 'First-Law Efficiency: Thermal vs. Electrical', 'FontWeight', 'bold');
legend(ax1, {'Thermal (\eta_{th})', 'Electrical (\eta_{elec})'}, 'Location', 'northwest');
for i = 1:4
    text(ax1, i, eta_tot_all(i) + 0.02 * max(eta_tot_all), sprintf('%.1f%%', eta_tot_all(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
end

ax2 = subplot(1, 3, 2);
set(ax2, 'Color', [1, 1, 1]);
b2 = bar(ax2, time_all, 'FaceColor', 'flat', 'BarWidth', 0.48);
b2.CData = bar_cols;
format_case_axis(ax2, top_lines, max(time_all) * 1.18);
ylabel(ax2, 'Warm-Up Time (Hours)', 'FontWeight', 'bold');
title(ax2, 'Time to Heat 200 L Tank (10 °C to 45 °C)', 'FontWeight', 'bold');
for i = 1:4
    text(ax2, i, time_all(i) + 0.02 * max(time_all), sprintf('%.1f h', time_all(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
end

ax3 = subplot(1, 3, 3);
set(ax3, 'Color', [1, 1, 1]);
b3 = bar(ax3, qloss_all, 'FaceColor', 'flat', 'BarWidth', 0.48);
b3.CData = bar_cols;
format_case_axis(ax3, top_lines, max(qloss_all) * 1.18);
ylabel(ax3, 'Mean Heat Loss (W)', 'FontWeight', 'bold');
title(ax3, 'Mean Radiative and Convective Heat Loss', 'FontWeight', 'bold');
for i = 1:4
    text(ax3, i, qloss_all(i) + 0.02 * max(qloss_all), sprintf('%.0f W', qloss_all(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
end

fig_traj = figure('Color', [1, 1, 1], 'InvertHardcopy', 'off', 'Position', [100, 100, 960, 720]);

ax_sub1 = subplot(2, 1, 1);
set(ax_sub1, 'Color', [1, 1, 1]);
hold(ax_sub1, 'on');
plot(ax_sub1, t_u_vec, qloss_u, 'LineWidth', 2.2, 'Color', [0.18, 0.63, 0.33]);
plot(ax_sub1, traj_pvt.time_h, traj_pvt.qloss, 'LineWidth', 2.2, 'Color', [0.12, 0.47, 0.71]);
plot(ax_sub1, traj_th.time_h,  traj_th.qloss,  'LineWidth', 2.2, 'Color', [0.84, 0.24, 0.24]);
plot(ax_sub1, traj_sel.time_h, traj_sel.qloss, 'LineWidth', 2.2, 'Color', [0.58, 0.20, 0.65]);
xlim(ax_sub1, [0, 20]);
if ~isempty(t_u_vec)
    idx_pv1 = find(t_u_vec <= 20, 1, 'last');
    if isempty(idx_pv1)
        idx_pv1 = 1;
    end
    text(ax_sub1, t_u_vec(idx_pv1), qloss_u(idx_pv1), '40 hrs total', ...
        'Color', [0.18, 0.63, 0.33], 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
end
hold(ax_sub1, 'off');
grid(ax_sub1, 'off');
box(ax_sub1, 'off');
xlabel(ax_sub1, 'Time (Hours)', 'FontWeight', 'bold');
ylabel(ax_sub1, 'Convective and Radiative Losses (W)', 'FontWeight', 'bold');
title(ax_sub1, 'Convective and Radiative Losses vs. Time Across Configurations', 'FontWeight', 'bold');
legend(ax_sub1, {'Unglazed PV', 'Glazed PVT', 'Pure Thermal (\epsilon_p = 0.9)', 'Pure Thermal (\epsilon_p = 0.3)'}, 'Location', 'best');

ax_sub2 = subplot(2, 1, 2);
set(ax_sub2, 'Color', [1, 1, 1]);
hold(ax_sub2, 'on');
plot(ax_sub2, t_u_vec, tp_u_vec, 'LineWidth', 2.2, 'Color', [0.18, 0.63, 0.33]);
plot(ax_sub2, traj_pvt.time_h, traj_pvt.Tp_c, 'LineWidth', 2.2, 'Color', [0.12, 0.47, 0.71]);
plot(ax_sub2, traj_th.time_h,  traj_th.Tp_c,  'LineWidth', 2.2, 'Color', [0.84, 0.24, 0.24]);
plot(ax_sub2, traj_sel.time_h, traj_sel.Tp_c, 'LineWidth', 2.2, 'Color', [0.58, 0.20, 0.65]);
xlim(ax_sub2, [0, 20]);
if ~isempty(t_u_vec)
    idx_pv2 = find(t_u_vec <= 20, 1, 'last');
    if isempty(idx_pv2)
        idx_pv2 = 1;
    end
    text(ax_sub2, t_u_vec(idx_pv2), tp_u_vec(idx_pv2), '40 hrs total', ...
        'Color', [0.18, 0.63, 0.33], 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
end
hold(ax_sub2, 'off');
grid(ax_sub2, 'off');
box(ax_sub2, 'off');
xlabel(ax_sub2, 'Time (Hours)', 'FontWeight', 'bold');
ylabel(ax_sub2, 'Plate Temperature, T_p (°C)', 'FontWeight', 'bold');
title(ax_sub2, 'Plate Temperature (T_p) vs. Time Across Configurations', 'FontWeight', 'bold');
legend(ax_sub2, {'Unglazed PV', 'Glazed PVT', 'Pure Thermal (\epsilon_p = 0.9)', 'Pure Thermal (\epsilon_p = 0.3)'}, 'Location', 'best');

function format_case_axis(ax, top_lines, y_top)
    xticks(ax, 1:4);
    set(ax, 'XTickLabel', {}, 'FontSize', 9.5);
    xlim(ax, [0.35, 4.65]);
    ylim(ax, [0, y_top]);
    grid(ax, 'off');
    box(ax, 'off');

    yl = ylim(ax);
    yspan = yl(2) - yl(1);
    for i = 1:4
        text(ax, i, yl(1) - 0.03 * yspan, top_lines{i}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'Interpreter', 'none', 'FontSize', 8.8, 'Clipping', 'off');
    end
end

function [m, traj] = run_glazed_metrics(p, M_tank, target_temp, dt, max_steps, P_pump, options)
    T_tank = 10 + 273.15;

    eth_h = zeros(max_steps, 1);
    eel_h = zeros(max_steps, 1);
    ql_h  = zeros(max_steps, 1);
    tp_h  = zeros(max_steps, 1);
    n_used = 0;

    for step = 1:max_steps
        p.Tin = T_tank;
        Tp_s = fzero(@(Tp) glazed_residual(Tp, p), T_tank + 10, options);
        [~, ~, Q_u_s, ~, Q_e_s, eta_e, Q_r_s, Q_c_s] = glazed_residual(Tp_s, p);

        n_used = step;
        eth_h(step) = (Q_u_s / (p.G * p.Ac)) * 100;
        eel_h(step) = eta_e * 100;
        ql_h(step)  = Q_r_s + Q_c_s;
        tp_h(step)  = Tp_s - 273.15;

        if T_tank >= target_temp
            break;
        end
        T_tank = T_tank + ((Q_u_s + (Q_e_s - P_pump)) * dt) / (M_tank * p.Cp);
    end

    m.eta_th = mean(eth_h(1:n_used));
    m.eta_el = mean(eel_h(1:n_used));
    m.qloss  = mean(ql_h(1:n_used));
    m.time_h = (n_used - 1) / 60;

    traj.time_h   = ((0:(n_used - 1))') / 60;
    traj.eta_comb = eth_h(1:n_used) + eel_h(1:n_used);
    traj.qloss    = ql_h(1:n_used);
    traj.Tp_c     = tp_h(1:n_used);
end

function [F1, Tout, Q_useful, Q_elec, Q_loss_conv, Q_loss_rad] = bare_residual(Tp, params)
    G = params.G;
    Tamb = params.Tamb;
    Tsky = params.Tsky;
    h_wind = params.h_wind;
    Tin = params.Tin;
    Ac = params.Ac;
    alpha = params.alpha;
    eps_p = params.eps_p;
    sigma = 5.670374419e-8;

    T_ref = 25 + 273.15;
    eta_ref = 0.25;
    beta_tc = 0.0026;

    eta_elec = max(0, eta_ref * (1 - beta_tc * (Tp - T_ref)));
    Q_elec = eta_elec * G * Ac;

    Tout = Tin;
    Q_useful = 0;

    Q_loss_conv = h_wind * Ac * (Tp - Tamb);
    Q_loss_rad  = eps_p * sigma * Ac * (Tp^4 - Tsky^4);

    Net_Solar_Thermal = (alpha * G * Ac) - Q_elec;
    F1 = Net_Solar_Thermal - Q_useful - (Q_loss_conv + Q_loss_rad);
end

function [F1, Tout, Q_useful, h_gap, Q_elec, eta_elec, Q_loss_rad, Q_loss_conv, Tg] = glazed_residual(Tp, params)
    G = params.G;
    Tamb = params.Tamb;
    Tsky = params.Tsky;
    h_wind = params.h_wind;
    Tin = params.Tin;
    mdot = params.mdot;
    Cp = params.Cp;
    Ac = params.Ac;
    alpha = params.alpha;
    eps_p = params.eps_p;
    P = params.P;
    L = params.L;
    L_gap = params.Lgap;
    hfi = params.hfi;
    sigma = 5.670374419e-8;

    err = 1;
    beta_tilt = deg2rad(45);
    Tg = (3 * Tp + h_wind * Tamb) / (3 + h_wind);

    while err > 1e-2
        Tg_old = Tg;
        T_mean = (Tp + Tg) / 2;
        k_air = 0.0241 * (T_mean / 273.15)^0.9;
        nu_air = 1.48e-5 * (T_mean / 293.15)^1.5;
        alpha_air = 2.88e-5 * (T_mean / 293.15)^1.75;
        beta_vol = 1 / T_mean;
        Ra = (9.81 * beta_vol * abs(Tp - Tg) * L_gap^3) / (nu_air * alpha_air);
        Ra_cos = max(Ra * cos(beta_tilt), 1e-9);
        term1_a = max(0, 1 - 1708 / Ra_cos);
        term1_b = 1 - (1708 * (sin(1.8 * beta_tilt))^1.6) / Ra_cos;
        term2   = max(0, (Ra_cos / 5830)^(1/3) - 1);
        Nu = 1 + 1.44 * term1_a * term1_b + term2;
        h_gap = (Nu * k_air) / L_gap;
        Tg = (h_gap * Tp + h_wind * Tamb) / (h_gap + h_wind);
        err = abs(Tg - Tg_old);
    end

    Fsky = (1 + cos(beta_tilt)) / 2;
    Fground = 1 - Fsky;
    T_ref = 298.15;
    eta_ref = 0.25;
    beta_tc = 0.0026;

    if isfield(params, 'run_elec') && ~params.run_elec
        eta_elec = 0;
        Q_elec = 0;
    else
        eta_elec = max(0, eta_ref * (1 - beta_tc * (Tp - T_ref)));
        Q_elec = eta_elec * G * Ac;
    end

    Tout = Tp - (Tp - Tin) * exp(-(P * L * hfi) / (mdot * Cp));
    Q_useful = mdot * Cp * (Tout - Tin);
    Q_loss_rad  = eps_p * sigma * Ac * (Fsky * (Tp^4 - Tsky^4) + Fground * (Tp^4 - Tamb^4));
    Q_loss_conv = h_gap * Ac * (Tp - Tg);
    Net_Solar_Thermal = (alpha * G * Ac) - Q_elec;
    F1 = Net_Solar_Thermal - Q_useful - (Q_loss_rad + Q_loss_conv);
end

function hfi = calc_hfi(mdot, D_i, mu, k_water, Cp)
    Re_D = 4 * mdot / (3.14 * D_i * mu);
    Pr = Cp * mu / k_water;
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
    hfi = k_water * Nu_D / D_i;
end