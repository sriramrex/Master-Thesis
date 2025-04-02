clear; close all; clc

% Given discrete natural frequencies (in Hz)
f_n_values = [150.89, 159.44, 165.59, 171.39, 171.39, 182.84, 189.69];

% Other system parameters
zeta = 0.1;      % Damping ratio

% Define excitation frequency range (in Hz)
f = 50:0.5:250;  

% Prepare the figure and color order
figure;
hold on;
colors = lines(length(f_n_values));  % Generates distinguishable colors

% Loop over each discrete natural frequency
for i = 1:length(f_n_values)
    f_n = f_n_values(i);
    
    % Calculate the isolator frequency for this mode
    f_i = f_n / sqrt(2);
    
    % Frequency ratio based on isolator frequency
    r = f / f_i;
    
    % Calculate displacement transmissibility:
    T_x = 1 ./ sqrt((1 - r.^2).^2 + (2*zeta*r).^2);

    % Max Transmissibility
    T_max = max(T_x);

    % Calculate the crossing point
    T_c = 1/sqrt((1 - (f_n/f_i)^2)^2 + (2*zeta*(f_n/f_i))^2);

    
    % Plot the transmissibility curve
    plot(f, T_x, 'LineWidth', 2, 'Color', colors(i,:));
    
    % Mark the isolator and natural frequencies on the plot
    plot(f_i, T_max, 'go','MarkerSize',10);
    plot(f_n, T_c, 'ro','MarkerSize',10);
end

% Enhance the plot appearance
xlabel('Excitation Frequency (Hz)');
ylabel('Displacement Transmissibility, T_x');
title('Displacement Transmissibility for Discrete Natural Frequencies');
grid on;

% Create legend entries for clarity
legendEntries = cell(1, length(f_n_values));
for i = 1:length(f_n_values)
    legendEntries{i} = sprintf('f_n = %.2f Hz', f_n_values(i));
end
legend(legendEntries, 'Location', 'Best');

hold off;
