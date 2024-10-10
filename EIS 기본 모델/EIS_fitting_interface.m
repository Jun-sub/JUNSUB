clc, clear, close all;

% assign data_path
path_data = 'C:\Users\admin\Documents\GitHub\JunSub\EIS 기본 모델\data_1.xlsx';
addpath 'C:\Users\admin\Documents\GitHub\JS_fitting_intro'
% load data
  data_e = importdata(path_data);

  f_e = data_e.data(:,2);         %[Hz]
  z_e_real = data_e.data(:,3);    %[Ohm]
  z_e_imag = data_e.data(:,4);    %[Ohm]
  z_e = z_e_real - 1i*z_e_imag;

%pre-plot
  figure(1)
  plot(z_e_real, -z_e_imag)
  legend('experimental');
% fitting

  % initial guess
  R_init = 1.56;
  C_init = 10^-3;
  A_init = 3;
  R2_init = 18.5;

  para_init = [R_init, C_init, A_init, R2_init];

% calculate initial guess
  Z_init = Z_model_RCW(f_e, para_init); %%f 값에 2pi 곱하지 않아도 되는건지 궁금합니다
%pre-plot2
  figure(1); 
  hold on 
  plot(real(Z_init), -imag(Z_init))
  legend('experimental', 'guess');

%Main fitting

%Cost function define
objfunc = @(para)func_cost(z_e,f_e,para);

%minimization
para_hat = fmincon(objfunc,para_init,[],[],[],[],[0,0,0,0],10*para_init);

%plot
Z_hat = Z_model_RCW(f_e,para_hat);

% pre-plot2
figure(1); hold on
plot(real(Z_hat), -imag(Z_hat))
xlabel('Re(Z)/Ohm');
ylabel('-Im(Z)/Ohm');
title('Impedance'); 
axis equal;
grid on;
legend('experimental', 'guess','fitting');

function cost = cost_func(z_e,f_e,para)
  z_mod = Z_model_RCW(f_e,para);
  cost = sum(sqrt((real(z_e - z_mod).^2))) + sum(sqrt((imag(z_e - z_mod).^2)));
%%sqrt 써야 하는지
end 
  