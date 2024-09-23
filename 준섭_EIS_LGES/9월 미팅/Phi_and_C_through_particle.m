clc; clear; close all;

x = linspace(0, 2, 1000); % x-axis representing the thickness direction (normalized 0 to 1)
frequency = 5; % Constant frequency for the sine wave
% amplitude_c = (1/3)*(1 - 0.3*x); % () determine the amplitude (Linear)

amplitude_c = zeros(size(x));

% amplitude_c = (3/2) * exp(-6 * x); % High frequency, freq 20 (Exp)
amplitude_c = (0.4) * exp(-0.5 * x); %low freq, freq = 5


% Create the concentration profile using the sine function
c_profile = amplitude_c .* sin(2 * pi * frequency * x);
% Plotting
figure;
plot(x, c_profile, 'LineWidth', 1.5,'Color', 'blue');hold on;
viscircles([2,0],2,'Color','black');
viscircles([5,0],1,'Color','black');
yline(0,'-k','LineWidth',2,'Alpha',0.3);
xlim([0,5])
ylim([-2,2])
xlabel('Particle thickness');
ylabel('Normalized value');
title('AC C_s profile at low freq');
% grid on;
box on;
legend('C_s', 'Active material particle')
hold off;