clc; clear; close all;

% Parameters
x = linspace(0, 1, 1000); % x-axis representing the thickness direction (normalized 0 to 1)
frequency = 25; % Constant frequency for the sine wave

% Define amplitude as a function of x, symmetrically with respect to x = 0.5
amplitude_c = zeros(size(x));

% Define amplitude based on conditions
% amplitude_c(x <= 0.3) = (0.8) * exp(-30 * (x(x <= 0.3))); %High freq,
% freq 80
% amplitude_c(x > 0.3 & x < 0.7) = (1/3) * exp(-20 * (x(x > 0.3 & x < 0.7)));
% amplitude_c(x >= 0.7) = (0.8) * exp(-30 * (1 - x(x >= 0.7)));

 amplitude_c(x <= 0.5) = (0.4) * exp(-6 * (x(x <= 0.5))); %low freq, freq 25
% amplitude_c(x > 0.3 & x < 0.7) = (1/3) * exp(-5 * (x(x > 0.3 & x < 0.7)));
 amplitude_c(x > 0.5) = (0.4) * exp(-6 * (1 - x(x >= 0.5)));



amplitude_p = (1/30);  % Constant amplitude for phi_profile

% Create the concentration profile using the sine function
c_profile = amplitude_c .* sin(2 * pi * frequency * x); 
phi_profile = amplitude_p .* sin(2 * pi * frequency * x);

% Plotting
figure;
plot(x, c_profile, 'LineWidth', 1.5, 'Color', 'blue'); hold on;
plot(x, phi_profile, 'LineWidth', 2, 'Color', 'red');

% Annotations for clarification
text(0.05, max(c_profile) + 0.2, 'Cathode', 'HorizontalAlignment', 'left', 'FontSize', 12);
text(0.95, max(c_profile) + 0.2, 'Anode', 'HorizontalAlignment', 'right', 'FontSize', 12);
text(0.5, min(c_profile) - 0.2, 'Separator', 'HorizontalAlignment', 'center', 'FontSize', 12);

% Labels and Title
xlabel('Electrode Thickness');
ylabel('Normalized value');
title('AC C_l & phi_l profile at low freq');
box on;
legend('C_l' , 'Phi_l');
ylim([-2,2]);

hold off;
