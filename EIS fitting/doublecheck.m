

% assign data_path
path_data = "C:\Users\gmdms\Documents\GitHub\Test\JS_EIS_practice\example1_eis_fc_soc48_t25.txt";

% load data
  data_e = importdata(path_data);

  f_e = data_e.data(:,1);         %[Hz]
  z_e_real = data_e.data(:,2);    %[Ohm]
  z_e_imag = -data_e.data(:,3);   %[Ohm]
  z_e = z_e_real + 1i*z_e_imag;

%pre-plot
  figure(1)
  plot(z_e_real, -z_e_imag)
  legend('experimental');
% fitting

  % initial guess
  R_init = 0.02;
  C_init = 1;
  A_init = 0.002;
  R2_init = 0.027;

  para_init = [R_init, C_init, A_init, R2_init];

% calculate initial guess
  Z_init = Z_model_RCW(f_e, para_init);

%pre-plot2
  figure(1); 
  hold on 
  plot(real(Z_init), -imag(Z_init))
  legend('experimental', 'guess');

%Main fitting

%Cost function define
objfunc = @(para)cost_func(z_e,f_e,para);

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

function cost = cost_func(z_e,f,para)
  z_mod = Z_model_RCW(f,para);
  cost = sum(real(z_e - z_mod).^2) + sum(imag(z_e - z_mod).^2);

end 
  