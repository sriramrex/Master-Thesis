clear; close all; clc

% Given system parameters
f_n = 135;       % Original system natural frequency in Hz
m = 1.51;        % Mass (kg)
% Choose a damping ratio for the isolator
zeta = 0.1;      

% Determine minimum isolator frequency (guideline)
f_i_min = f_n / sqrt(2);  % ~95.5 Hz
fprintf('Minimum isolator frequency (guideline): %.2f Hz\n', f_i_min);

% For illustration, let's design an isolator at this frequency
f_i = f_i_min;  
omega_i = 2*pi*f_i;  % Convert to rad/s if needed

% Define excitation frequency range (in Hz)
f = 50:0.5:200;  
r = f / f_i;       % Frequency ratio

% Calculate displacement transmissibility
T_x = 1 ./ sqrt((1 - r.^2).^2 + (2*zeta*r).^2);

% Max Transmissibility

T_max = max(T_x);

% Calculate the crossing point

T_c = 1/sqrt((1 - (f_n/f_i)^2)^2 + (2*zeta*(f_n/f_i))^2);

% Plot the transmissibility curve
figure;
plot(f, T_x, 'LineWidth',2);
xlabel('Excitation Frequency (Hz)');
ylabel('Displacement Transmissibility, T_x');
xline(f_i_min, '--b', 'Modal Frequency', 'LineWidth', 1.5)
title('Displacement Transmissibility vs Excitation Frequency');
grid on;

% Mark the operating point at f = 135 Hz
hold on;
plot(f_i_min, T_max, 'bo', 'MarkerSize',10);
plot(f_n, T_c, 'ro', 'MarkerSize',10);
xline(f_n, '--r', 'Natural Frequency', 'LineWidth', 1.5)
legend('T_x','Operating Frequency (135 Hz)');
hold off;
