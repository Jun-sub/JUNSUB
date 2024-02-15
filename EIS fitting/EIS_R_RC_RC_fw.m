clear all; clc; close;

f_min = -3;
f_max = 5;
points = 1000;

f = logspace(f_min, f_max, points);
w = 2* pi* f;

R0 = 0.5;

R1 = 1;
R2 = 0.5;

C1 = 10^-5;
C2 = 10^-3;
A = 2;
t = 1; % t = time constant (tau)




%Calculation
Z_W = A .* coth(sqrt(1i *w .*t))./sqrt(1i *w .*t); %Warburg impedance
Z_C1 = 1 ./(1i*w*C1);
Z_C2 = 1 ./(1i*w*C2);

%Total impedacne
Z = R0 + (R1 .*Z_C1) ./(R1 + Z_C1) + (R2 .*Z_C2) ./(R2 + Z_C2) + Z_W;

%plot
figure;
plot(real(Z), -imag(Z));
xlabel('Re(Z)/Ohm');
ylabel('-Im(Z)/Ohm');
xlim([0,3]);
ylim([0,3]);
title('Impedacne');
grid on
