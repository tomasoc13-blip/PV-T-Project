clc, close all, clear

% Time (hours)
time = 0:23;

% Global irradiance (W/m^2)
G0 = [0 0 0 0 0 0 10 76 161 250 325 372 ...
      394 381 339 269 181 92 19 0 0 0 0 0];

G20 = [0 0 0 0 0 0 9 83 189 302 398 457 ...
       489 472 417 327 214 103 17 0 0 0 0 0];

G45 = [0 0 0 0 0 0 7 84 201 329 439 506 ...
       544 525 462 358 230 106 13 0 0 0 0 0];

G70 = [0 0 0 0 0 0 5 74 186 308 415 480 ...
       517 499 438 337 213 95 8 0 0 0 0 0];

% Smooth curve interpolation
time_fine = linspace(0,23,500);

G0_fit  = interp1(time,G0,time_fine,'pchip');
G20_fit = interp1(time,G20,time_fine,'pchip');
G45_fit = interp1(time,G45,time_fine,'pchip');
G70_fit = interp1(time,G70,time_fine,'pchip');

% Plot
figure;

plot(time_fine,G0_fit,'r-','LineWidth',2, ...
    'DisplayName','0°');
hold on;

plot(time_fine,G20_fit,'k-','LineWidth',2, ...
    'DisplayName','20°');

plot(time_fine,G45_fit,'b-','LineWidth',2, ...
    'DisplayName','45°');

plot(time_fine,G70_fit,'m-','LineWidth',2, ...
    'DisplayName','70°');

% Original data points
plot(time,G0,'rs','MarkerSize',5,'HandleVisibility','off');
plot(time,G20,'kd','MarkerSize',5,'HandleVisibility','off');
plot(time,G45,'bo','MarkerSize',5,'HandleVisibility','off');
plot(time,G70,'m^','MarkerSize',5,'HandleVisibility','off');

% Axis formatting
xlabel('Time (hrs)');
ylabel('Global Irradiance (W/m^2)');

xticks(0:3:23);
xticklabels({'00:45','03:45','06:45','09:45', ...
    '12:45','15:45','18:45','21:45'});

xlim([0 23]);
title('PVGIS: Aran Islands Daily Irradiance in September (Varying Tilt-Angle)')

legend('Location','northwest','FontWeight','bold');
set(gca,'FontSize',14,'FontWeight','bold');
grid off;
box off;

% 45 degree data: average within +/- 4 hours of peak

% Find peak
[G_peak, idx] = max(G45);
t_peak = time(idx);

% Integration limits
t_start = t_peak - 4;
t_end = t_peak + 4;

% Fine time vector within integration region
t_avg = linspace(t_start,t_end,1000);

% Interpolate irradiance
G_avg_curve = interp1(time,G45,t_avg,'pchip');

% Numerical integration
G_integral = trapz(t_avg,G_avg_curve);

% Average irradiance
G_average = G_integral/(t_end-t_start);

fprintf('Peak irradiance = %.1f W/m^2\n',G_peak);
fprintf('Peak time = %.2f hours\n',t_peak);
fprintf('Integration period = %.2f to %.2f hours\n',t_start,t_end);
fprintf('Average irradiance = %.2f W/m^2\n',G_average);