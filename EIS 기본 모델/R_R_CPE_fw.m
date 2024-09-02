clear all; , clc; close;

f_min = -3;
f_max = 3;
points = 1000;

f = logspace(f_min, f_max, points);
w = 2* pi* f;

R = 1;
C = 10^-3;
A = 3;
t = 200; % t = time constant (tau)
Q = 1;
a = 5; % a stands for alpha



%Calculation
Z_W = A .* coth(sqrt(1i *w .*t))./sqrt(1i *w .*t); %Warburg impedance
Z_CPE = 1 ./(Q .*(1i .*w).^a);

%Total impedacne
Z = (R.*Z_CPE)./(R + Z_CPE) + Z_W;

%plot
figure;
plot(real(Z), -imag(Z));
xlim([0,2]);
ylim([0,2]);
xlabel('Re(Z)/Ohm');
ylabel('-Im(Z)/Ohm');
title('Impedacne');
grid on
