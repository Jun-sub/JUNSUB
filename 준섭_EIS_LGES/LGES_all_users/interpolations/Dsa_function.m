function [ Ds_est ] = Dsa_function(theta0a,T)
% This function interpolates the solid-state diffusivity from a measured data
% Data from "Samsung22V\n_am1_Ds_HC1067.txt"

% Currently only for T = 25 degC;

data = [0.954433145	2.71235E-15
        0.507974605	2.14359E-14
        0.123300993	2.55E-13];

Ds0 = interp1(data(:,1),data(:,2),theta0a,'linear','extrap');
Ea = 30e3; % [J/mol] = 30 KJ/mol
R = 8.314; % [J/mol.K]
Ds_est = Ds0*exp(-Ea/R*(1/T - 1/298.15));

end

