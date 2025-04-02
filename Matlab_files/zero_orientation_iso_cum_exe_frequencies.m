clear; close all; clc

% Discrete excitation frequencies (Hz)
f_excitation = [150.89, 159.44, 165.59, 171.39, 171.39, 182.84, 189.69];

% Damping ratio for the isolators
zeta = 0.1;

% Define a broad range of frequencies for plotting the transmissibility curves
f = 50:0.5:250;  

% Create a figure and set up a color map for distinct curves
figure;
hold on;
colors = lines(length(f_excitation));

% Preallocate an array for curve handles
curveHandles = gobjects(length(f_excitation),1);

for i = 1:length(f_excitation)
    % Current excitation frequency
    f_exc = f_excitation(i);
    % Calculate the corresponding isolator frequency
    f_iso = f_exc / sqrt(2);
    
    % Compute frequency ratio based on isolator frequency
    r = f / f_iso;
    
    % Calculate displacement transmissibility:
    T_x = 1 ./ sqrt((1 - r.^2).^2 + (2*zeta*r).^2);

    % Max Transmissibility
    T_max = max(T_x);
    
    % Plot the transmissibility curve and store its handle
    curveHandles(i) = plot(f, T_x, 'LineWidth', 2, 'Color', colors(i,:));
    
    % Mark the excitation frequency point on the curve (exclude from legend)
    r_exc = f_exc / f_iso;
    T_x_exc = 1 / sqrt((1 - r_exc^2)^2 + (2*zeta*r_exc)^2);
    plot(f_exc, T_x_exc, 'o', 'Color', colors(i,:), 'MarkerSize', 8, 'MarkerFaceColor', colors(i,:), 'HandleVisibility', 'off');
    plot(f_iso, T_max, 'o','MarkerSize',10,'LineWidth', 2);
    % Plot a vertical dashed line at the isolator frequency (exclude from legend)
    xline(f_iso, '--', 'Color', colors(i,:), 'LineWidth', 1, 'HandleVisibility', 'off');
end

% Enhance the plot with labels and a title
xlabel('Excitation Frequency (Hz)');
ylabel('Displacement Transmissibility, T_x');
title('Transmissibility Curves for Given Excitation Frequencies and Their Isolator Frequencies');
grid on;

% Create a legend for the transmissibility curves only
legendEntries = cell(1, length(f_excitation));
for i = 1:length(f_excitation)
    legendEntries{i} = sprintf('f_{exc} = %.2f Hz, f_{iso} = %.2f Hz', f_excitation(i), f_excitation(i)/sqrt(2));
end
legend(curveHandles, legendEntries, 'Location', 'Best');

hold off;
