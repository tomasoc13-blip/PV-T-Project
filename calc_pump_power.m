function [P_pump, dP_total, Re_D] = calc_pump_power(V_litre_min, L, D_i)

% Inputs:
%   V_litre_min : Volumetric flow rate (L/min)
%   L           : Pipe length (m)
%   D_i         : Pipe inner diameter (m)
%   eta_pump    : pump efficiency 
%
% Outputs:
%   P_pump      : Required electrical pump power (W)
%   dP_total    : Total pressure drop including 1.5 FoS for bends/fittings (Pa)
%   Re_D        : Reynolds number 

    
eta_pump = 0.70;


    % Guard against zero or negative flow
    if V_litre_min <= 0
        P_pump = 0;
        dP_total = 0;
        Re_D = 0;
        return;
    end

    % Fluid properties (Water at ~25 °C, matching your existing model)
    ro = 998;                  % Density (kg/m^3)
    mu = 0.890e-3;             % Dynamic viscosity (Pa*s)
    FoS = 1.5;                 % Factor of safety for bends, elbows, and fittings

    % Convert flow rate and compute flow velocity
    V_dot = V_litre_min / 60000;        % m^3/s
    A_pipe = (pi / 4) * D_i^2;          % m^2
    v_flow = V_dot / A_pipe;            % m/s

    %  Reynolds number
    Re_D = (ro * v_flow * D_i) / mu;

   
 % Surface roughness for drawn copper tubing (m)
eps_roughness = 1.5e-6; % 0.0015 mm

% 3. Darcy friction factor (f)
if Re_D < 2300
    f_darcy = 64 / Re_D; % Laminar flow (independent of roughness)
else
    % Haaland  correlation for turbulent flow in rough pipes
    rel_roughness = eps_roughness / D_i;
    inv_sqrt_f = -1.8 * log10( (rel_roughness / 3.7)^1.11 + (6.9 / Re_D) );
    f_darcy = 1 / (inv_sqrt_f^2);
end

    % Major pressure drop + 1.5 FoS for minor losses (bends/fittings)
    dP_major = f_darcy * (L / D_i) * (0.5 * ro * v_flow^2);
    dP_total = FoS * dP_major;          % Pa

    % Required electrical pump power
    P_pump = (dP_total * V_dot) / eta_pump; % Watts
end