clear; close all; clc

% Define system parameters
f_n = 90;              % Natural frequency in Hz
omega_n = f_n * 2 * pi; % Convert natural frequency to rad/s
m = 1.51;               % Mass (kg)
zeta = 0.1;             % Damping ratio (0.01 = 10% damping)

% Frequency range obtained duing the manual testing of the tool
f = 1:0.5:300;         % Excitation frequency in Hz
omega = f * 2 * pi;     % Convert to rad/s

% Calculate stiffness
k = m * omega_n^2;

% Frequency ratio (r = omega/omega_n)
r = omega ./ omega_n;


%% Transmissibility Calculations

% Force Transmissibility 
t_f = sqrt((1 + (2*zeta*r).^2) ./ sqrt((1 - r.^2).^2 + (2*zeta*r).^2));

% Displacement Amplitude Ratio (X/X0 - for forced vibration)

% X/X0 = 1 / sqrt((1 - r.^2).^2 + (2*zeta*r).^2)
t_x = 1 ./ sqrt((1 - r.^2).^2 + (2*zeta*r).^2);

% Acceleration Transmissibility

t_a = (r.^2)./sqrt((1 - r.^2).^2 + (2*zeta*r).^2);


%% Plotting

figure('Position', [100, 100, 1200, 500])

% Plot Force Transmissibility
subplot(1,2,1)
plot(f, t_f, 'b', 'LineWidth', 2)
hold on
xline(f_n, '--r', 'Natural Frequency', 'LineWidth', 1.5);
xlabel('Excitation Frequency (Hz)')
ylabel('Force Transmissibility (TR)')
title('Force Transmissibility vs Frequency')
grid on
%axis tight
legend('Transmissibility', 'Location', 'northeast')

% Plot Displacement Ratio
subplot(1,2,2)
plot(f, t_x, 'k', 'LineWidth', 2)
hold on
xline(f_n, '--r', 'Natural Frequency', 'LineWidth', 1.5);
xlabel('Excitation Frequency (Hz)')
ylabel('Displacement Ratio (X/X_0)')
title('Displacement Amplitude Ratio vs Frequency')
grid on
axis tight
legend('Displacement Ratio', 'Location', 'northeast')

% Plot Acceleration Ratio
% subplot(1,3,3)
% semilogy(f, t_a, 'r', 'LineWidth', 2)
% hold on
% xline(f_n, '--r', 'Natural Frequency', 'LineWidth', 1.5);
% xlabel('Excitation Frequency (Hz)')
% ylabel('Acceleration Ratio')
% title('Acceleration Amplitude Ratio vs Frequency')
% grid on
% axis tight
% legend('Acceleration Ratio', 'Location', 'northeast')




% ==================================================================
% Additional Notes
% ==================================================================
% For steady-state displacement calculations:
% F0 = ...;       % Define your input force magnitude
% phase_angle = atan2(2*zeta*r, 1 - r.^2);
% x_steady = (F0/k) * displacement_ratio .* sin(omega.*t - phase_angle);



