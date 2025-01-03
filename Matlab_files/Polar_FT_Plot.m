clear; close all; clc

base_path = ['/home/sduddu-iit.local/HHCM_LAB/Matlab_Files/Isolated_Tool_Test/FTSensor0_streamingData_09_12_2024_12_19_03.csv'];
ftsensor = readtable(base_path);
sensor_t = ftsensor.SampleID;
sensor_data = [ftsensor.Fx'; ftsensor.Fy'];

%% Compute the polar plot for the x,y forces

time_stamp = linspace(0,15000,length(sensor_t));
time = time_stamp/1000;
period = 0.001; % Time for one complete cycle (2*pi)
% Calculate angles
time_angle = 2 * pi * mod(time, period) / period;

% force values 

f_x = ftsensor.Fx';
f_y = ftsensor.Fy';



% Calculate  euclidien distance
distance = sqrt(f_x.^2 + f_y.^2);

% Calculate angles
angles = atan2(f_y,f_x);
%angles = wrapTo360(theta);


%% 3d data

x = distance .*cos(angles);
y = distance .*sin(angles);

x_t = distance .*cos(time_angle);
y_t = distance .*sin(time_angle);

%% Plot distance vs angle
figure;
polarscatter(angles, distance,'filled'); % Polar plot for distance as a function of angle
title('Distance vs Angle');

figure(2);

%plot3(f_x(1:500),f_y(1:500),time(1:500),"filled");
scatter3(x,y,time,36,time,"filled");
colormap jet;
colorbar;
% Label the axes
xlabel('X (Radius * cos(Angle))');
ylabel('Y (Radius * sin(Angle))');
zlabel('Time');
title('3D Polar Plot with Radius, Angle, and Time');

figure(3);

polarscatter(time_angle, distance,"filled"); % Polar plot for distance as a function of angle
title('Distance vs Angle');

figure(4);

scatter3(x_t,y_t,time,36,time,"filled");
colormap jet;
colorbar;
% Label the axes
xlabel('X_time ');
ylabel('Y_time ');
zlabel('Time');
title('3D Polar Plot with Radius, Angle, and Time')

figure (5);

subplot(3,1,1);
plot(time, distance,'b','LineWidth',2);
xlabel('time');
ylabel('radius');

subplot(3,1,2);
plot(time, angles,'g','LineWidth',2);
xlabel('time');
ylabel('angle (rad)');

subplot(3,1,3);
plot(time,time_angle,'r','LineWidth',2);
xlabel('time');
ylabel('time @ angle (rad)');







