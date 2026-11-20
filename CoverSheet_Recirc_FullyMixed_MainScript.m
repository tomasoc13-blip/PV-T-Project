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
