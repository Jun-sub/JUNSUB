
clear; clc; close all;

%% Configuration

%EIS data path
    folder = '';
    file = '';

% SOC and T
    soc = 0.5; %[1]
    T = 298.15; %[K]

%Fitting configuration
    type_weight = 1; % 0 for absolute error, 1 for relative error
    type_acf = 2; % 1 for anode, 2 for cathode, 3 for full cell

% Optimization options
    options = optimset('display', 'iter', 'MaxIter', 100, 'MaxFunEvals', 1e5, ...
        'Tolfun',1e-8,'TolX',1e-8,'FinDifftype','central');

    % iter: 반복마다 최적화 정보 표시 설정
    % MaxIter, 100: 최대 반복 횟수 설정 (100번)
    % MaxFunEvals, 1e5: 함수 평가의 최대 횟수 설정
    % Tolfun, 1e-8: 함수 값의 변화에 대한 허용 오차, 1e-8 이하가 되면 수렴으로 간주.
    % TolX, 1e-8: 최적화 변수의 변화에 대한 허용 오차, 1e-8 이하가 되면 수렴으로 간주.
    % FinDifftype, cnetral: 유한 차분 방식(Finite difference method)를 사용해 미분을 계산.
    % central은 중앙 차분 방식 (central different method) 


%Parameters
    bounds = [
        0.5 2 %(1) R_itsc
        0.1 50; % (2) i0
        0.1 10; % (3) C_dl
        0.01 10; % (4) Ds
        0.1 10; %(5) kappa_el
        0.001 10; % (6) D_el
        0.1 100; % (7) av
        ];

    lb = bounds(:,1); % lower bounds
    up = bounds(:,2); % upper bounds

    factors_ini = [1 1 1 1 1 1 1]; %Factor*bounds 식에서 bounds 변화를 통해서 계산.

%% Load and Pre-processing Data

    % load EIS data
    data = csvread([folder filesep file]);
    f_data = data(:,1);
    z_re_data = data(:,2);
    z_im_data = data(:,3); 
    z_data = [z_re_data z_im_data];

    figure(1) %pre-plot
    plot(z_re_data, -z_im_data, 'o'); 
    hold on;
    grid on;

    axis_limit = 1.1*max(max(abs(z_data))); %데이터 절대값 --> 열방향 최대 값 두 개 --> 행방향 최대값
    set(gca,'Box', 'on',... % Box: 축을 네 모서리에 모두 표시 
    'PlotBoxAspectRatio',[1 1 1],... % 플롯의 가로, 세로, 깊이 비율 설정. [1 1 1]의 경우는 정사각형(정육면체) 플롯 박스의 비율을 설정
    'FontUnits', 'Points', 'FontSize', 10, 'FontName', 'Times New Roman', ... %FontUnits, Points: 폰트 단위를 포인트로 설정; FontSize, 10: 폰트 크기를 10 포인트로 설정; FontName, Times New Roman: 폰트 이름을 Times New Roman으로 설정
    'XLim', [0 axis_limit], 'YLim', [0 axis_limit]); %x축, y축 범위 설정. 
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')
    
    % trim high-frequency inductance part 
    f_data = f_data(z_im_data <= 0);
    z_re_data = z_re_data(z_im_data <= 0);
    z_im_data = z_im_data(z_im_data <= 0);
    z_data = [z_re_data z_im_data];
    figure(1)
    plot(z_re_data, -z_im_data,'o');
    hold on;
    grid on;

    %% Define weighting vector
        if type_weight == 1 % Minimizing the relative error
            weight = (z_re_data.^2 + z_im_data.^2).^(-0.5); %제곱 후 역수 --> 큰 수의 가중치를 작게, 작은 수의 가중치를 크게
            weight_matrix = [weight weight];
        elseif type_weight == 0 % Minimizinf the absolute error
            weight_matrix = ones(size(z_re_data));
        end


    %% Plot Results
  %  [z_model0, paras_used0] = BSL_func_EISmodel(f_data,factors_ini, soc, T, type_acf)
  %  [z_model1, paras_used1] = BSL_func_EISmodel(f_data,factors_hat, soc, T, type_acf) 
  % 이 부분 아직 이해 불가. 

    %% Call EIS model
   %  weighted_model = @(factors,f_data)BSL_func_EISmodel(f_data, factors,soc,T,type_acf)...
   %     .*weight_matrix;
   % weighted_data = z_data.*weight_matrix;
   % 
   % tic;
   % [factors_hat, resnorm,residual,~,~,~,jacobian_hat] ...
   %      = lsqcurvefit(weighted_model,factors_ini,...
   %                    f_data,weighted_data, lb, ub, options);
   % toc;


    % Nyquist plot
    figure(2)
    plot(z_data(:,1), -z_data(:,2), 'ok', 'linewidth', 1);
    hold on;
    plot(z_model0(:,1), -z_model0(:,2), 'ob', 'linewidth', 1);
    % plot(z_model1(:,1), -z_model1(:,2), 'ob', 'linewidth', 1);
    legend('Exp data', 'Model Initail'); % 'Model Fit'
    daspect ([1 1 2]) %데이터 자체의 축 비율을 설정

    axis_limit = 1.1*max(max(abs(z_data)));
    set(gca,'Box','on',... %위와 동일한 설정. 
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... 
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    hold off
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')

    % Zoom-in semicircle
    figure(3) %우선 그냥 플랏
    plot(z_data(:,1), -z_data(:,2), 'ok', 'linewidht', 1); hold on
    plot(z_model0(:,1), -z_model0(:,2), 'ob', 'linewidht', 1);
    plot(z_model1(:,1), -z_model1(:,2), 'or', 'linewidht', 1);
    legend('Exp Data', 'Model Initial', 'Model fit')

    f_zoom_lb = 10; %[Hz]
    idx_zoom = f_data > f_zoom_lb %10 Hz보다 큰 데이터의 인덱스 설정
    axis_limit = 1.1*max(max(abs(z_data(idx_zoom))));
    set(gca,'Box','on',... %위와 동일한 설정. 
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... 
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    hold off
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')

    %% Result Summary

    Result.factors_hat = factors_hat';
    Result.paras_hat = paras_used1';
    Result.z_model = z_model1; 