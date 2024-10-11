clc, clear, close all;

%Frequency
w = logspace(-3,3,300);

% initial guess (3RC-W)
R0 = 3;

C1 = 10^-3;
R1 = 10;
A1 = 0;

C2 = 10^-2;
R2 = 20;
A2 = 0;

C3 = 10^-1;
R3 = 30;
A3 = 0;


%% Randle
%Warburg impedance
Z_W1 = A1 .* (1-1i) ./ sqrt(w);
Z_W2 = A2 .* (1-1i) ./ sqrt(w);
Z_W3 = A3 .* (1-1i) ./ sqrt(w);

%Impedance components
Z_RW1 = R1 + Z_W1;
Z_C1 = 1 ./ (1i*w*C1);

Z_RW2 = R2 + Z_W2;
Z_C2 = 1 ./ (1i*w*C2);

Z_RW3 = R3 + Z_W3;
Z_C3 = 1 ./ (1i*w*C3);

%Total impedance
Z1 = R0 + (Z_RW1 .* Z_C1) ./ (Z_RW1 + Z_C1) + (Z_RW2 .* Z_C2) ./ (Z_RW2 + Z_C2) + (Z_RW3 .* Z_C3) ./ (Z_RW3 + Z_C3);

%Plot
figure(1); hold on;
plot(real(Z1), -imag(Z1),'ro-','LineWidth',1.5)

%% RC + Zw
%Warburg impedance
Z_W1 = A1 .* (1-1i) ./ sqrt(w);
Z_W2 = A2 .* (1-1i) ./ sqrt(w);
Z_W3 = A3 .* (1-1i) ./ sqrt(w);

%Impedance components
Z_R1 = R1;
Z_C1 = 1 ./ (1i*w*C1);

Z_R2 = R2;
Z_C2 = 1 ./ (1i*w*C2);

Z_R3 = R3;
Z_C3 = 1 ./ (1i*w*C3);

%Total impedance
Z2 = R0 + (1./Z_R1 + 1./Z_C1).^-1 + (1./Z_R2 + 1./Z_C2).^-1 + (1./Z_R3 + 1./Z_C3).^-1 + (1./Z_W1 + 1./Z_W2 + 1./Z_W3).^-1;

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
