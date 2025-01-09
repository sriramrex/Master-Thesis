clear; close all; clc

base_path = 'C:\Users\sduddu\SensLiveView\FTSensor0_streamingData_03_01_2025_14_47_27.csv';
ftsensor = readtable(base_path);
sensor_t = ftsensor.SampleID;
sensor_data = [ftsensor.Fx'; ftsensor.Fy'];

%% Compute the polar plot for the x,y forces

time_stamp = linspace(0,20000,length(sensor_t));
time = time_stamp/1000;
period = 20; % Time for one complete cycle (2*pi)
% Calculate angles
time_angle = 2 * pi * mod(time, period) / period;

% force values 

Fx = ftsensor.Fx';
Fy = ftsensor.Fy';
Fz = ftsensor.Fz';

% Calculate polar coordinates (magnitude and angle)
force_magnitude = sqrt(Fx.^2 + Fy.^2);
%force_angle = atan2(Fy, Fx); % Angle in radians

force_angle = mod(atan2(Fy, Fx), 2*pi);


% Plot in polar coordinates
figure;
polarscatter(force_angle, force_magnitude);
title('Polar Distribution of Fx and Fy');

% Additional visualization (histogram of angles)
figure;
histogram(rad2deg(force_angle), 50);
title('Angle Distribution of Fx and Fy');
xlabel('Angle (degrees)');
ylabel('Frequency of Angles');

