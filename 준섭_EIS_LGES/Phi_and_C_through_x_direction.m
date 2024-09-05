clc; clear; close all;

% Parameters
x = linspace(0, 1, 1000); % x-axis representing the thickness direction (normalized 0 to 1)
frequency = 60; % Constant frequency for the sine wave

% Define amplitude as a function of x (larger near the electrodes, smaller near the separator)
amplitude_c = (1/8)*(1 - x/0.5); % () determine the amplitude
amplitude_p = (1/30);
% Create the concentration profile using the sine function
c_profile = amplitude_c .* sin(2 * pi * frequency * x);
phi_profile = amplitude_p .* sin(2 * pi * frequency * x);
% Plotting
figure;
plot(x, c_profile, 'LineWidth', 2,'Color','blue');hold on;
plot(x, phi_profile, 'LineWidth', 2,'Color','red');


% Annotations for clarification
text(0.05, max(c_profile) + 0.2, 'Cathode', 'HorizontalAlignment', 'left', 'FontSize', 12);
text(0.95, max(c_profile) + 0.2, 'Anode', 'HorizontalAlignment', 'right', 'FontSize', 12);
text(0.5, min(c_profile) - 0.2, 'Separator', 'HorizontalAlignment', 'center', 'FontSize', 12);

% Labels and Title
xlabel('Electrode Thickness');
ylabel('Normalized value');
title('AC C_l & phi_l profile at high freq');
% grid on;
box on;
legend('C_l' , 'Phi_l')
ylim([-2,2])
hold off;