clear; clc; close all;

%% Configurations

% EIS data path
path_folder = 'C:\Users\admin\Desktop\켄텍 서류\과제 폴더\LGES\데이터\modeling_DC and EIS\12_6cm2_soc10_EIS # Sample 2';
name_file = 'PEIS_C11_cathode_cycle_soc30.csv';

% Load EIS data
data = csvread([path_folder filesep name_file]);
f_data = data(:, 1);
z_re_data = data(:, 2);
z_im_data = data(:, 3);
z_data = [z_re_data z_im_data];

figure(1);
plot(z_re_data, -z_im_data, 'o'); hold on; grid on;

%% Warburg include
% Parameters for the RC circuit and Warburg impedance
R0 = 0.08; % Setting R0 to zero to start from origin
R = 0.2;  % Resistance in ohms
C = 1e-4; % Capacitance in farads

% Frequency range (logarithmically spaced)
f = logspace(-3, 6, 100); % Frequency from 0.01 Hz to 1 MHz
omega = 2 * pi * f; % Angular frequency

% Calculate the impedance of the RC circuit
Z_RC_re = R0 + R ./ (1 + (omega .* R .* C).^2); % Real part of RC
Z_RC_im = -(omega .* R.^2 .* C) ./ (1 + (omega .* R .* C).^2); % Imaginary part of RC

% Total impedance
Z_total_re = Z_RC_re;
Z_total_im = Z_RC_im;

% Plot the Nyquist plot
figure(1);
plot(Z_total_re, -Z_total_im, 'ko'); % Plotting the semi-circle with Warburg
hold on;
xlabel("Z'");
ylabel("Z''");
grid on;
axis equal;
xlim([0 0.3]);
ylim([0 0.3]);
legend('Experimental', 'Initial guess');
hold off;