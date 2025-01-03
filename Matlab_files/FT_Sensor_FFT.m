clear; close all; clc

% % Define the base path and file names for the 7 FT sensor readings
% base_path = '/home/sduddu-iit.local/HHCM_LAB/Matlab_Files/Isolated_Tool_Test/';
% file_names = { ...
%     'FTSensor0_streamingData_09_12_2024_12_19_03.csv', ...
%     'FTSensor0_streamingData_09_12_2024_12_21_17.csv', ...
%     'FTSensor0_streamingData_09_12_2024_12_23_32.csv', ...
%     'FTSensor0_streamingData_09_12_2024_12_25_59.csv', ...
%     'FTSensor0_streamingData_09_12_2024_12_28_10.csv', ...
%     'FTSensor0_streamingData_09_12_2024_12_29_46.csv', ...
%     'FTSensor0_streamingData_09_12_2024_12_31_43.csv'};

% % %This is for the Free_Free_FT_Sensor readings
% base_path = '/home/sduddu-iit.local/HHCM_LAB/Matlab_Files/Free_Free_Cond_FTsensor/';
% file_names = { ...
%     'FTSensor0_streamingData_25_11_2024_13_02_12.csv', ...
%     'FTSensor0_streamingData_25_11_2024_14_20_14.csv', ...
%     'FTSensor0_streamingData_25_11_2024_14_28_19.csv', ...
%     'FTSensor0_streamingData_25_11_2024_14_31_54.csv', ...
%     'FTSensor0_streamingData_25_11_2024_14_35_18.csv', ...
%     'FTSensor0_streamingData_25_11_2024_14_43_50.csv', ...
%     'FTSensor0_streamingData_25_11_2024_14_45_30.csv'};

% % These are the reading of the Force Torque Sensor in O orientation.
% base_path = '/home/sduddu-iit.local/HHCM_LAB/Matlab_Files/FTSensor_Test_Results';
% file_names = { ...
%     'FTSensor0_streamingData_20_12_2024_15_39_30.csv', ...
%     'FTSensor0_streamingData_20_12_2024_15_41_19.csv', ...
%     'FTSensor0_streamingData_20_12_2024_15_52_11.csv', ...
%     'FTSensor0_streamingData_20_12_2024_16_00_42.csv', ...
%     'FTSensor0_streamingData_20_12_2024_16_15_48.csv', ...
%     'FTSensor0_streamingData_20_12_2024_16_25_51.csv', ...
%     'FTSensor0_streamingData_20_12_2024_16_32_30.csv'};



% % These are the reading of the Force Torque Sensor in 9O orientation.
% base_path = '/home/sduddu-iit.local/HHCM_LAB/Matlab_Files/FTSensor_Test_Results';
% file_names = { ...
%     'FTSensor0_streamingData_20_12_2024_17_27_34.csv', ...
%     'FTSensor0_streamingData_20_12_2024_18_04_21.csv', ...
%     'FTSensor0_streamingData_20_12_2024_17_50_19.csv', ...
%     'FTSensor0_streamingData_20_12_2024_17_58_44.csv', ...
%     'FTSensor0_streamingData_20_12_2024_18_13_48.csv', ...
%     'FTSensor0_streamingData_20_12_2024_18_18_27.csv', ...
%     'FTSensor0_streamingData_20_12_2024_18_32_42.csv'};


% % These are the reading of the Force Torque Sensor in 45 orientation.
% base_path = '/home/sduddu-iit.local/HHCM_LAB/Matlab_Files/FTSensor_Test_Results';
% file_names = { ...
%     'FTSensor0_streamingData_20_12_2024_12_21_22.csv', ...
%     'FTSensor0_streamingData_20_12_2024_12_55_20.csv', ...
%     'FTSensor0_streamingData_20_12_2024_12_58_03.csv', ...
%     'FTSensor0_streamingData_20_12_2024_13_02_25.csv', ...
%     'FTSensor0_streamingData_20_12_2024_13_04_44.csv', ...
%     'FTSensor0_streamingData_20_12_2024_13_07_27.csv', ...
%     'FTSensor0_streamingData_20_12_2024_13_13_25.csv'};


% These are the reading of the Force Torque Sensor in 9O orientation (Updated).
base_path = '/home/sduddu-iit.local/HHCM_LAB/Matlab_Files/FTSensor_Test_Results';
file_names = { ...
    'FTSensor0_streamingData_23_12_2024_18_29_01.csv', ...
    'FTSensor0_streamingData_23_12_2024_18_32_49.csv', ...
    'FTSensor0_streamingData_20_12_2024_17_50_19.csv', ...
    'FTSensor0_streamingData_23_12_2024_18_37_20.csv', ...
    'FTSensor0_streamingData_23_12_2024_18_48_14.csv', ...
    'FTSensor0_streamingData_23_12_2024_18_49_38.csv', ...
    'FTSensor0_streamingData_23_12_2024_18_52_13.csv'};


filelength = length(file_names);

% Initialize an empty table to store results for all sensors
all_sensor_data = table();

% Mass of the System
mass = 1.51; % Mass of the tool in Kg.

% Arrays to store natural frequency and stiffness values for all sensors
natural_freq_x_all = zeros(filelength, 1);
natural_freq_y_all = zeros(filelength, 1);
natural_freq_z_all = zeros(filelength, 1);

stiffness_x_all = zeros(filelength, 1);
stiffness_y_all = zeros(filelength, 1);
stiffness_z_all = zeros(filelength, 1);

% High-pass filter parameters (e.g., cutoff frequency of 1 Hz)
cutoff_frequency = 100; % Set cutoff frequency in Hz
filter_order = 10; % Order of the filter

% Loop through all the FT sensor data files
for i = 1:filelength
    
    % Read the FT sensor data for the current file
    file_path = fullfile(base_path, file_names{i});
    FTSensor_1 = readtable(file_path);
    sensor_t_1 = FTSensor_1.SampleID;
    sensor_data_1 = [FTSensor_1.Fx'; FTSensor_1.Fy'; FTSensor_1.Fz'];

    %% Compute FFT for all the axes
    time_stamp = linspace(0, 20000, length(sensor_t_1));
    time = time_stamp / 1000; 
    dt = time(2) - time(1); % Time Difference
    Fs = 1 / dt;  % Sampling Frequency
    L = length(time); % Number of samples

    % High-pass filter design (using Butterworth filter)
    [b, a] = butter(filter_order, cutoff_frequency / (Fs / 2), 'high'); % High-pass filter coefficients

    % Filter the data for each axis
    filtered_data_x = filter(b, a, sensor_data_1(1,:));
    filtered_data_y = filter(b, a, sensor_data_1(2,:));
    filtered_data_z = filter(b, a, sensor_data_1(3,:));

    %% FFT on filtered data

    % FFT for X-axis Force (Filtered)
    f_x = fft(filtered_data_x);
    P2_x = abs(f_x / L);           % Two-sided spectrum
    P1_x = P2_x(1:L / 2 + 1);      % Single-sided spectrum
    P1_x(2:end - 1) = 2 * P1_x(2:end - 1); % Double non-DC components
    f = Fs * (0:(L / 2)) / L;      % Frequency range

    % FFT for Y-axis Force (Filtered)
    f_y = fft(filtered_data_y);
    P2_y = abs(f_y / L);
    P1_y = P2_y(1:L / 2 + 1);
    P1_y(2:end - 1) = 2 * P1_y(2:end - 1);

    % FFT for Z-axis Force (Filtered)
    f_z = fft(filtered_data_z);
    P2_z = abs(f_z / L);
    P1_z = P2_z(1:L / 2 + 1);
    P1_z(2:end - 1) = 2 * P1_z(2:end - 1);



    %% Dominant frequency (approximate natural frequency) for each axis

    thre_x = max(P1_x)*0.3;% Set threshold to 30% of the maximum amplitude
    idx_x = find(P1_x>thre_x, 1);
    dominant_frequency_x = f(idx_x);

    thre_y = max(P1_y)*0.3;
    idx_y = find(P1_y>thre_y, 1);
    dominant_frequency_y = f(idx_y); 

    thre_z = max(P1_z)*0.3;
    idx_z = find(P1_z>thre_z, 1);
    dominant_frequency_z = f(idx_z);


    % Dominant frequency (approximate natural frequency) for each axis
    % [~, idx_x] = max(P1_x); % Find the peak in the X-axis
    % dominant_frequency_x = f(idx_x); 
    % 
    % [~, idx_y] = max(P1_y); % Find the peak in the Y-axis
    % dominant_frequency_y = f(idx_y); 
    % 
    % [~, idx_z] = max(P1_z); % Find the peak in the Z-axis
    % dominant_frequency_z = f(idx_z);


    % Stiffness Calculation (N/m) for each axis
    omega_n_x = 2 * pi * dominant_frequency_x; % Natural frequency in rad/s for X-axis
    stiffness_x = mass * omega_n_x^2;  % Stiffness for X-axis

    omega_n_y = 2 * pi * dominant_frequency_y; % Natural frequency in rad/s for Y-axis
    stiffness_y = mass * omega_n_y^2;  % Stiffness for Y-axis

    omega_n_z = 2 * pi * dominant_frequency_z; % Natural frequency in rad/s for Z-axis
    stiffness_z = mass * omega_n_z^2;  % Stiffness for Z-axis

    % Store the stiffness values in arrays
    natural_freq_x_all(i) = dominant_frequency_x ;
    natural_freq_y_all(i) = dominant_frequency_y ;
    natural_freq_z_all(i) = dominant_frequency_z ;

    stiffness_x_all(i) = stiffness_x;
    stiffness_y_all(i) = stiffness_y;
    stiffness_z_all(i) = stiffness_z;

    %% Calculating the essential data such as RMS, min, max of the Forces

    % For X-axis
    rms_x = sqrt(mean(filtered_data_x.^2));
    mean_x = mean(filtered_data_x);
    max_x = max(filtered_data_x);
    min_x = min(filtered_data_x);

    % For Y-axis
    rms_y = sqrt(mean(filtered_data_y.^2));
    mean_y = mean(filtered_data_y);
    max_y = max(filtered_data_y);
    min_y = min(filtered_data_y);

    % For Z-axis
    rms_z = sqrt(mean(filtered_data_z.^2));
    mean_z = mean(filtered_data_z);
    max_z = max(filtered_data_z);
    min_z = min(filtered_data_z);


    %% Plots for each FT sensor data

    % % Create a new figure for each reading
    % figure;
    % 
    % % X-axis plot
    % subplot(3, 1, 1);
    % plot(f, P1_x);
    % hold on;
    % plot(dominant_frequency_x, P1_x(idx_x), 'ro', 'MarkerFaceColor', 'r'); % Mark natural frequency
    % text(dominant_frequency_x, P1_x(idx_x), sprintf('  %.2f Hz', dominant_frequency_x), 'Color', 'red');
    % title(['Amplitude Spectrum of Force (X-axis) - Sensor ', num2str(i)]);
    % xlabel('Frequency (Hz)');
    % ylabel('|Amplitude|');
    % grid on;
    % legend('Amplitude Spectrum', 'Natural Frequency');
    % 
    % % Y-axis plot
    % subplot(3, 1, 2);
    % plot(f, P1_y);
    % hold on;
    % plot(dominant_frequency_y, P1_y(idx_y), 'ro', 'MarkerFaceColor', 'r'); % Mark natural frequency
    % text(dominant_frequency_y, P1_y(idx_y), sprintf('  %.2f Hz', dominant_frequency_y), 'Color', 'red');
    % title(['Amplitude Spectrum of Force (Y-axis) - Sensor ', num2str(i)]);
    % xlabel('Frequency (Hz)');
    % ylabel('|Amplitude|');
    % grid on;
    % legend('Amplitude Spectrum', 'Natural Frequency');
    % 
    % % Z-axis plot
    % subplot(3, 1, 3);
    % plot(f, P1_z);
    % hold on;
    % plot(dominant_frequency_z, P1_z(idx_z), 'ro', 'MarkerFaceColor', 'r'); % Mark natural frequency
    % text(dominant_frequency_z, P1_z(idx_z), sprintf('  %.2f Hz', dominant_frequency_z), 'Color', 'red');
    % title(['Amplitude Spectrum of Force (Z-axis) - Sensor ', num2str(i)]);
    % xlabel('Frequency (Hz)');
    % ylabel('|Amplitude|');
    % grid on;
    % legend('Amplitude Spectrum', 'Natural Frequency');

    %% Tabular display for each sensor reading
    ForceData = table({'X-axis'; 'Y-axis'; 'Z-axis'}, ...
                      [rms_x; rms_y; rms_z], ...
                      [mean_x; mean_y; mean_z], ...
                      [max_x; max_y; max_z], ...
                      [min_x; min_y; min_z], ...
                      [dominant_frequency_x; dominant_frequency_y; dominant_frequency_z], ...
                      [stiffness_x; stiffness_y; stiffness_z], ...
                      'VariableNames', {'Axis', 'RMS (N)', 'Mean (N)', 'Max (N)', 'Min (N)', 'Natural Frequency (Hz)', 'Stiffness (N/m)'});

    % Display the final table in the command window
    disp(['Combined Results for Speed ' num2str(i) ':']);
    disp(ForceData);
   
    % Append the current sensor data to the all_sensor_data table
    all_sensor_data = [all_sensor_data; ForceData];
   
end

% Save the entire table to a CSV file
writetable(all_sensor_data, 'FTSensor_Analysis_Results.csv');

%% Visualization: Plot the stiffness variation across different sensors
figure;

% Grouped bar plot for stiffness values
bar_data_stiff = [stiffness_x_all, stiffness_y_all, stiffness_z_all];
bar(1:filelength, bar_data_stiff, 'grouped');
xlabel('Sensor Number');
ylabel('Stiffness (N/m)');
title('Stiffness Variation Across Different Speeds');
legend('X-axis', 'Y-axis', 'Z-axis','Location','northwest');
grid on;

figure;

% Grouped bar plot for natural frequency values
bar_data_freq = [natural_freq_x_all, natural_freq_y_all, natural_freq_z_all];
% Grouped bar plot for natural frequency values
b = bar(1:filelength, bar_data_freq, 'grouped'); % Create the grouped bar plot

% Set colors for each group
b(1).FaceColor = 'r'; % Red for X-axis
b(2).FaceColor = 'g'; % Green for Y-axis
b(3).FaceColor = 'b'; % Blue for Z-axis

% Add labels and title
xlabel('Tool Speeds');
ylabel('Natural Frequency (Hz)');
legend('X-axis', 'Y-axis', 'Z-axis','Location','northwest');
title('Natural Frequency Across Different Speeds');

% Annotate the bars with the natural frequency values
for i = 1:filelength
    % Annotate X-axis bar
    text(i - 0.22, bar_data_freq(i, 1) + 2, num2str(natural_freq_x_all(i), '%.2f'),'FontWeight' , 'Bold','Color', 'r' , 'HorizontalAlignment', 'center');
    
    % Annotate Y-axis bar
    text(i, bar_data_freq(i, 2) + 2, num2str(natural_freq_y_all(i), '%.2f'), 'FontWeight' , 'Bold','Color', 'g', 'HorizontalAlignment', 'center');
    
    % Annotate Z-axis bar
    text(i + 0.22, bar_data_freq(i, 3) + 2, num2str(natural_freq_z_all(i), '%.2f'),'FontWeight' , 'Bold', 'Color', 'b', 'HorizontalAlignment', 'center');
end

grid on;
