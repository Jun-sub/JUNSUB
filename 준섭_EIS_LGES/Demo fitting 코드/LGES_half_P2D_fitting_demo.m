%% 
% This function 
% (1) finds the best-fit parameters for a halfcells
% calls a EIS model function

clear; clc; close all


%% Configurations

% EIS data path
    path_folder = 'G:\공유 드라이브\BSL-Data\LGES\LG raw data\12_6cm2_soc10_EIS # Sample 1';
    %path_folder = 'G:\Shared drives\BSL-Data\LGES\LG raw data\12_6cm2_soc10_EIS # Sample 2';
    name_file = 'PEIS_C09_cathode_cycle_soc90.csv';

    addpath('C:\Users\admin\Documents\GitHub\BSL_LGES2024')
% SOC and T (for initial guess - they are functions of soc, T)
    soc = 0.9; % [1]
    T = 298.15; %[K]

% Fitting configuration
    type_weight = 1; % 0 for absolute error, 1 for relative error
    type_acf =2; % 1 for anode, 2 for cathode, 3 for full cell

% Optimization options
    %move into for iteration


% Parameters 
    bounds = [...
         0.5 2 % (1) R_itsc
         0.1 50; % (2) i0
         0.1 50; % (3) C_dl
         0.1 50; % (4) Ds
         0.1 10; % (5) kappa_el
         0.01 10; % (6) D_el
         0.1 10; % (7) Av
         ]; 
    lb = bounds(:,1);
    ub = bounds(:,2);
    %factors_ini = [1 1 1 1 1 1 1]
    factors_ini = [0.8652 1 5.7685 1 2.0947 0.9428 0.5346]; %only for io &
    %Ds

%% Load and Pre-processing Data

    % load EIS data
    data = load([path_folder filesep name_file]);
    f_data = data(:,1);
    z_re_data = data(:,2);
    z_im_data = data(:,3);
    z_data = [z_re_data z_im_data];

    figure(1)
    plot(z_re_data,-z_im_data,'o'); hold on; grid on
    
    axis_limit = 1.1*max(max(abs(z_data)));
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')

    % trim the high-frequency inductance part
    f_data = f_data(z_im_data<=0);
    z_re_data = z_re_data(z_im_data<=0);
    z_im_data = z_im_data(z_im_data<=0);
    z_data = [z_re_data z_im_data];
    figure(1)
    plot(z_re_data,-z_im_data,'o'); hold on; grid on


%% Define weighting vector
    
    if type_weight == 1 % if minimizing the relative error
    weight = (z_re_data.^2 + z_im_data.^2).^(-0.5);
    weight_matrix = [weight weight];    
    elseif type_weight == 0 % if minimizing the absolute error
    weight_matrix = ones(size(z_re_data));
    end




%% Call EIS model
factors_hat = factors_ini;


for n = 1:10 %53
options= optimset('display','Iter','MaxIter',1,'MaxFunEvals',1e5,...
        'TolFun',1e-8,'TolX',1e-8,'FinDiffType','central');

   weighted_model = @(factors,f_data)BSL_func_EISmodel_V1_half(f_data, factors,soc,T,type_acf)...
       .*weight_matrix;
   weighted_data = z_data.*weight_matrix;

   tic;
   [factors_hat, resnorm,residual,~,~,~,jacobian_hat] ...
        = lsqcurvefit(weighted_model,factors_hat,...
                      f_data,weighted_data, lb, ub, options);
   


%% Plot Results
[z_model0, paras0] = BSL_func_EISmodel_V1_half(f_data,factors_ini,soc,T,type_acf);
[z_model1, paras1] = BSL_func_EISmodel_V1_half(f_data,factors_hat,soc,T,type_acf);

% Nyquist Plot
figure(2)
plot(z_data(:,1),-z_data(:,2),'ok','linewidth',1); hold on
%plot(z_model0(:,1),-z_model0(:,2),'ob','linewidth',1)
plot(z_model1(:,1),-z_model1(:,2),'or','linewidth',1)
legend('Exp Data','P2D')
daspect ([1 1 2])


    axis_limit = 1.1*max(max(abs(z_data)));
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    set(gcf,'position',[200 700 600 600])
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')
    hold off
    


% Zoom-in semicircle
figure(3)
plot(z_data(:,1),-z_data(:,2),'ok','linewidth',1); hold on
%plot(z_model0(:,1),-z_model0(:,2),'ob','linewidth',1)
plot(z_model1(:,1),-z_model1(:,2),'or','linewidth',1)
legend('Exp Data','P2D')
daspect ([1 1 2])

    f_zoom_lb = 10; %[Hz] 
    idx_zoom = f_data>f_zoom_lb;
    axis_limit = 1.1*max(max(abs(z_data(idx_zoom,:))));
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    
    set(gcf,'position',[820 700 600 600])
    %pause(1)
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')
    hold off
  pause(2)
toc;
end 

fprintf ('Model was fitted to the experimental data')

% Opt_data = [paras0([2 4]') paras1([2 4]')] 
% openvar("Opt_data")
%% Result Summary
   
Result.factors_hat = factors_hat;
Result.paras_hat = paras1;
Result.z_model = z_model1;

