clear all; , clc; close;

f_min = -3;
f_max = 3;
points = 1000;

f = logspace(f_min, f_max, points);
w = 2* pi* f;

R = 1;
C = 10^-3;
A = 5;
t = 1; % t = time constant (tau)




%Calculation
Z_W = A .* coth(sqrt(1i *w .*t))./sqrt(1i *w .*t); %Warburg impedance
Z_RW = R + Z_W; % Real impedance + W
Z_C = 1 ./(1i*w*C);

%Total impedacne
Z = (Z_RW.*Z_C)./(Z_RW + Z_C);

%plot
figure;
plot(real(Z), -imag(Z));
xlim([0,3]);
ylim([0,3]);
xlabel('Re(Z)/Ohm');
ylabel('-Im(Z)/Ohm');
title('Impedacne');
grid on
