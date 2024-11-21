clc, clear, close all;

%Frequency
w = logspace(-3,3,300);

% initial guess
R0 = 3;
C = 10^-2;
A = 1.2;
R1 = 35;

%% Randle
%Warburg impedance
Z_W = A .* (1-1i) ./ sqrt(w);

%Impedance components
Z_RW = R1 + Z_W;
Z_C = 1 ./ (1i*w*C);

%Total impedance
Z1 = R0 + (Z_RW .* Z_C) ./ (Z_RW + Z_C);

%Plot
figure(1); hold on;
plot(real(Z1), -imag(Z1),'ro-','LineWidth',1.5)

%% RC + Zw
%Warburg impedance
Z_W = A .* (1-1i) ./ sqrt(w);

%Impedance components
Z_R = R1;
Z_C = 1 ./ (1i*w*C);

%Total impedance
Z2 = R0 + (1./Z_R + 1./Z_C).^-1 + Z_W;

%Plot
plot(real(Z2), -imag(Z2),'bo-','LineWidth',1.5)
hold off;

xmax = real(max(Z1));
xlim([0 1.1*xmax]);
ylim([0 1.1*xmax]);
grid on;
box on;
legend('Randle','RC + W')

cost = norm(Z1-Z2);
disp(cost);