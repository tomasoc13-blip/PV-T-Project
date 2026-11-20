function [F1, Tout, Q_useful, h_gap, Q_elec, eta_elec, Q_loss_rad, Q_loss_conv, Tg] = glazed_residual(Tp, params)
    % parameters
    G = params.G; Tamb = params.Tamb; Tsky = params.Tsky;
    h_wind = params.h_wind; Tin = params.Tin; mdot = params.mdot;
    Cp = params.Cp; Ac = params.Ac; alpha = params.alpha;
    eps_p = params.eps_p; P = params.P; L = params.L; L_gap = params.Lgap;
    hfi = params.hfi; 
    sigma = 5.670374419e-8; 
    err = 1;

    %  DYNAMIC h_gap CALCULATION (Hollands et al.) 
    beta_tilt = deg2rad(45); 

% Initial guess for Tg
Tg = (3 * Tp + h_wind * Tamb) / (3 + h_wind); % initial guess, h_gap = 3

while err > 1e-2   % Converge Tg 
    Tg_old = Tg;
    
    T_mean = (Tp + Tg) / 2;
    k_air = 0.0241 * (T_mean / 273.15)^0.9;
    nu_air = 1.48e-5 * (T_mean / 293.15)^1.5;
    alpha_air = 2.88e-5 * (T_mean / 293.15)^1.75;
    beta_vol = 1 / T_mean;
    
    Ra = (9.81 * beta_vol * abs(Tp - Tg) * L_gap^3) / (nu_air * alpha_air);
    Ra_cos = Ra * cos(beta_tilt);
    term1_a = max(0, 1 - 1708 / Ra_cos);
    term1_b = 1 - (1708 * (sin(1.8 * beta_tilt))^1.6) / Ra_cos;
    term2   = max(0, (Ra_cos / 5830)^(1/3) - 1);
    Nu = 1 + 1.44 * term1_a * term1_b + term2;
    
    h_gap = (Nu * k_air) / L_gap;
    Tg = (h_gap * Tp + h_wind * Tamb) / (h_gap + h_wind); % Series Resistance to compute Tg
    
    err = abs(Tg - Tg_old);
end

    Fsky = (1+cos(beta_tilt))/2; % Sky view factor, account for tilt angle
    Fground = 1 - Fsky;
    
    % PV Electrical Efficiency & Power
    T_ref = 298.15; 
    eta_ref = 0.25; beta_tc = 0.0026;

  
    eta_elec = eta_ref * (1 - beta_tc * (Tp - T_ref));
   
    
    Q_elec = eta_elec * G * Ac;
    %Q_elec = 0; % Pure thermal model

    
    Tout = Tp - (Tp - Tin) * exp(-(P * L * hfi) / (mdot * Cp));
    Q_useful = mdot * Cp * (Tout - Tin);
    
    Q_loss_rad  = eps_p * sigma * Ac *(Fsky * (Tp^4 - Tsky^4) + Fground*(Tp^4 - Tamb^4));
    Q_loss_conv = h_gap * Ac * (Tp - Tg);
    
    Net_Solar_Thermal = (alpha * G * Ac) - Q_elec;
    F1 = Net_Solar_Thermal - Q_useful - (Q_loss_rad + Q_loss_conv);
end