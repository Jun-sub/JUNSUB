%% 
% This function 
% (1) finds the best-fit parameters for a halfcells
% calls a EIS model function

clear; clc; close all

%% Configurations

% EIS data path
    path_folder = 'G:\공유 드라이브\BSL-Data\LGES\LG raw data\12_6cm2_soc10_EIS # Sample 1';
    name_file = 'PEIS_C09_cathode_cycle_soc90.csv';

    addpath('C:\Users\admin\Documents\GitHub\BSL_LGES2024')

% SOC and T (for initial guess - they are functions of soc, T)
    soc = 0.9; % [1]
    T = 298.15; %[K]

% Fitting configuration
    type_weight = 1; % 0 for absolute error, 1 for relative error
    type_acf = 2; % 1 for anode, 2 for cathode, 3 for full cell

% Optimization options
    %move into for iteration

% Parameters 
    bounds = [...
         0.5 2;   % (1) R_itsc
         0.1 50;  % (2) i0
         0.1 50;  % (3) C_dl
         0.1 1.2; % (4) Ds
         0.1 10;  % (5) kappa_el
         0.01 10; % (6) D_el
         0.1 10;  % (7) Av
         ]; 
    lb = bounds(:,1);
    ub = bounds(:,2);
    
    % Initial guess for optimization
    factors_ini = [0.9194 1 1.3086 1 1.9113 1.5774 4.2363]; %only for i0 & Ds
% Surf resolution
point = 50; %much point will take longer time
%% Load and Pre-processing Data

    % Load EIS data
    data = load([path_folder filesep name_file]);
    f_data = data(:,1);
    z_re_data = data(:,2);
    z_im_data = data(:,3);
    z_data = [z_re_data z_im_data];

    % Plot experimental data
    figure(1)
    plot(z_re_data, -z_im_data, 'o'); hold on; grid on
    
    axis_limit = 1.1 * max(max(abs(z_data)));
    set(gca, 'Box', 'on', ...
              'PlotBoxAspectRatio', [1 1 1], ...
              'FontUnits', 'points', 'FontSize', 10, 'FontName', 'Times New Roman', ...
              'XLim', [0 axis_limit], 'Ylim', [0 axis_limit])
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')

    % Trim the high-frequency inductance part
    f_data = f_data(z_im_data <= 0);
    z_re_data = z_re_data(z_im_data <= 0);
    z_im_data = z_im_data(z_im_data <= 0);
    z_data = [z_re_data z_im_data];
    plot(z_re_data, -z_im_data, 'o'); hold on; grid on

%% Define weighting vector
    if type_weight == 1 % if minimizing the relative error
        weight = (z_re_data.^2 + z_im_data.^2).^(-0.5);
        weight_matrix = [weight weight];    
    elseif type_weight == 0 % if minimizing the absolute error
        weight_matrix = ones(size(z_re_data));
    end

%% Call EIS model
factors_hat = factors_ini;

% Prepare for optimization tracking
i0_range = linspace(0.1, 1, point);
Ds_range = linspace(0.1, 1, point);
[i0_grid, Ds_grid] = meshgrid(i0_range, Ds_range);
resnorm_grid = zeros(size(i0_grid)); % Store resnorm values

% Initial calculation for resnorm grid
weighted_model = @(factors, f_data) BSL_func_EISmodel_V1_half(f_data, factors, soc, T, type_acf) ...
    .* weight_matrix;

% Calculate resnorm for each i0 and Ds value
for i = 1:length(i0_range)
    for j = 1:length(Ds_range)
        test_factors = [factors_ini(1), i0_range(i), factors_ini(3), Ds_range(j), ...
                        factors_ini(5), factors_ini(6), factors_ini(7)]; % Fix other parameters
        [test_model, paras_res] = BSL_func_EISmodel_V1_half(f_data, test_factors, soc, T, type_acf);
        resnorm_grid(j, i) = norm(z_data-test_model); % Calculate residual norm
    end
end

% resnorm_path = zeros(1,7);
% factor_path = zeros(7,2);


input('Press enter to start')

for n = 1:7 % Number of iterations
    options = optimset('display', 'iter', 'MaxIter', 1, 'MaxFunEvals', 1e5, ...
                       'TolFun', 1e-8, 'TolX', 1e-8, 'FinDiffType', 'central');

    weighted_data = z_data .* weight_matrix;

    tic;
    [factors_hat, resnorm, residual, ~, ~, ~, jacobian_hat] ...
        = lsqcurvefit(weighted_model, factors_hat, f_data, weighted_data, lb, ub, options);

    
    %fartor_hat trim
    factors_hat = [factors_ini(1), factors_hat(2), factors_ini(3), factors_hat(4), ...
                        factors_ini(5), factors_ini(6), factors_ini(7)];
    %% Plot Results
    [z_model0, paras0] = BSL_func_EISmodel_V1_half(f_data, factors_ini, soc, T, type_acf);
    [z_model1, paras1] = BSL_func_EISmodel_V1_half(f_data, factors_hat, soc, T, type_acf);

    % Nyquist Plot
    figure(2)
    plot(z_data(:,1), -z_data(:,2), 'ok', 'linewidth', 1); hold on
    plot(z_model1(:,1), -z_model1(:,2), 'or', 'linewidth', 1)
    legend('Exp Data', 'P2D')
    daspect([1 1 2])

    axis_limit = 1.1 * max(max(abs(z_data)));
    set(gca, 'Box', 'on', ...
              'PlotBoxAspectRatio', [1 1 1], ...
              'FontUnits', 'points', 'FontSize', 10, 'FontName', 'Times New Roman', ...
              'XLim', [0 axis_limit], 'Ylim', [0 axis_limit])
    set(gcf, 'position', [200 700 600 600])
    box on;
    grid on;
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')
    hold off

    % Zoom-in semicircle
    figure(3)
    plot(z_data(:,1), -z_data(:,2), 'ok', 'linewidth', 1); hold on
    plot(z_model1(:,1), -z_model1(:,2), 'or', 'linewidth', 1)
    legend('Exp Data', 'P2D')
    daspect([1 1 2])

    f_zoom_lb = 10; %[Hz] 
    idx_zoom = f_data > f_zoom_lb;
    axis_limit = 1.1 * max(max(abs(z_data(idx_zoom,:))));
    set(gca, 'Box', 'on', ...
              'PlotBoxAspectRatio', [1 1 1], ...
              'FontUnits', 'points', 'FontSize', 10, 'FontName', 'Times New Roman', ...
              'XLim', [0 axis_limit], 'Ylim', [0 axis_limit])
    box on;
    grid on;
    set(gcf, 'position', [820 700 600 600])
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')
    hold off
    
    % Resnorm Calculation
    resnorm_now = norm(z_data - z_model1);
    resnorm_path(n,1) = resnorm_now;
    factor_path(n,:) = [factors_hat(2),factors_hat(4)];
    % Resnorm visualization
    
    figure(4);
    cla;
    hold on;

    surf(i0_grid, Ds_grid, resnorm_grid);
    xlabel('i0');
    ylabel('Ds');
    zlabel('resnorm');
    title('Optimization Process Visualization');
    colorbar;
    view(120,50)
    box on;
    grid on;
    plot3(factor_path(:,1),factor_path(:,2),resnorm_path(:,1), 'ro', 'MarkerSize', 5, 'LineWidth',2,'MarkerFaceColor','red');
    plot3(factor_path(:,1),factor_path(:,2),resnorm_path(:,1), 'r-', 'LineWidth',2,'MarkerFaceColor','none');
    plot3(factor_path(n,1),factor_path(n,2),resnorm_path(n), 'co', 'MarkerSize', 6, 'LineWidth',2,'MarkerFaceColor','cyan');
    set(gcf, 'position', [1440 700 600 600])
    hold off;
    % plot3(i0_grid(min_i0_idx, min_Ds_idx), Ds_grid(min_i0_idx, min_Ds_idx), min_resnorm, 'ro', 'MarkerSize', 10, 'LineWidth', 3);

    pause(2)
    toc;
end 


fprintf('Model was fitted to the experimental data\n');



%% Result Summary
Result.factors_hat = factors_hat;
Result.paras_hat = paras1;
Result.z_model = z_model1;
