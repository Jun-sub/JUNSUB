function [ ] = plotting(  )
%PLOTTING Summary of this function goes here
%   Detailed explanation goes here


N = 100;
soc_vec = linspace(0,1,N);

xa_max = 0.8781;                            xc_max = 0.9319;                              % {modified}
xa_min = 0.0216;                            xc_min = 0.3532;                              % {modified}
xa_vec  = xa_max-(1-soc_vec)*(xa_max-xa_min);  xc_vec = xc_min+(1-soc_vec)*(xc_max-xc_min);     % {modified}

%Ua_chg_vec = Uc_function_v2(xa_vec,1);
% Uc_dchg_vec = Uc_function_v2(xa_vec,0);

dUcdx_dchg_vec = -(Uc_function_v2(xc_vec+0.01,0) -  Uc_function_v2(xc_vec-0.01,0))/0.002;
dUadx_dchg_vec = -(Ua_function_v2(xa_vec+0.01,0) -  Ua_function_v2(xa_vec-0.01,0))/0.002;

figure(1)
subplot(2,1,1); plot(soc_vec,dUcdx_dchg_vec)
subplot(2,1,2); plot(soc_vec,dUadx_dchg_vec)

end

