function [F1, Tout, Q_useful, Q_elec, Q_loss_conv, Q_loss_rad] = bare_residual(Tp, params)
    % Unpack environmental parameters
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
    hfi = params.hfi;
    sigma = 5.670374419e-8; 

    % --- PV Electrical Efficiency Calculations ---
    T_ref = 25 + 273.15;      % Reference temperature (25 °C in Kelvin)
    eta_ref = 0.25;           % Base electrical efficiency (e.g., 25%)
    beta_tc = 0.0026;         % Temperature coefficient
    
    % Actual efficiency drops as plate temperature (Tp) rises
    eta_elec = eta_ref * (1 - beta_tc * (Tp - T_ref));
    eta_elec = max(0, eta_elec); 
    
    % Electrical power generated (W)
    Q_elec = eta_elec * G * Ac;
    
    % 1. Calculate T_out using the constant wall temperature HX equation
    Tout = Tp - (Tp - Tin) * exp(-(P * L * hfi) / (mdot * Cp));
    
    % 2. Calculate useful thermal heat removed by water
    Q_useful = mdot * Cp * (Tout - Tin);
    
    % 3. Calculate top losses (convection to wind + radiation to sky)
    Q_loss_conv = h_wind * Ac * (Tp - Tamb);
    Q_loss_rad  = eps_p * sigma * Ac * (Tp^4 - Tsky^4);
    
    % 4. The Residual
    Net_Solar_Thermal = (alpha * G * Ac) - Q_elec;
    F1 = Net_Solar_Thermal - Q_useful - (Q_loss_conv + Q_loss_rad);
end