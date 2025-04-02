clear; close all; clc

% Define system parameters
w_n = 135;        % Natural frequency (Hz or rad/s as appropriate)
m   = 1.51;       % Mass (kg)
zeta = 0.1;       % Damping ratio

% Stiffness calculation
k = m * w_n^2;

% Frequency range for analysis
w = 85:0.5:300;   % Testing frequencies

% Normalize frequency (r = ω/ωₙ)
r = w / w_n;

% Calculate displacement transmissibility (ratio)
T_x = 1 ./ sqrt((1 - r.^2).^2 + (2*zeta*r).^2);

% Calculate force transmissibility (if base excitation is applied)
T_f = sqrt((1 + (2*zeta*r).^2) ./ ((1 - r.^2).^2 + (2*zeta*r).^2));

% Calculate acceleration transmissibility
T_a = (r.^2) ./ sqrt((1 - r.^2).^2 + (2*zeta*r).^2);

% For steady state displacement plot:
% Assuming a unit base excitation amplitude A = 1.
A = 1;
x_p = A * T_x;  % Steady state response amplitude

% Plotting all four graphs in one figure with subplots
figure;

% 1. Displacement Transmissibility
subplot(2,2,1);
plot(w, T_x, 'LineWidth',2);
xlabel('Frequency');
ylabel('Displacement Transmissibility');
title('Displacement Transmissibility vs Frequency');
grid on;

% 2. Force Transmissibility
subplot(2,2,2);
plot(w, T_f, 'LineWidth',2);
xlabel('Frequency');
ylabel('Force Transmissibility');
title('Force Transmissibility vs Frequency');
grid on;

% 3. Acceleration Transmissibility
subplot(2,2,3);
plot(w, T_a, 'LineWidth',2);
xlabel('Frequency');
ylabel('Acceleration Transmissibility');
title('Acceleration Transmissibility vs Frequency');
grid on;

% 4. Steady State Displacement
subplot(2,2,4);
plot(w, x_p, 'LineWidth',2);
xlabel('Frequency');
ylabel('Steady State Displacement');
title('Steady State Displacement vs Frequency');
grid on;
