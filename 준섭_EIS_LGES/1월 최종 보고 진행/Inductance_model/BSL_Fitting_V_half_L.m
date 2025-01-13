%% 
% This function 
% (1) finds the best-fit parameters for a halfcells calls a EIS model function
% (2) improve fitting result by addapting additional distribution model

%Configuration 설정 후 최적화 진행

%결과는 para_sum으로 저장, (1) SOC P2D model 사용 (2) R_itsc (3) i0 (4) C_dl (5) Ds (6) kappa_el (7) D_el (8) Av
%P2D + Dist model 사용 (11) R_itsc (12) i0 (13) C_dl (14) Ds (15) kappa_el
%(16) D_el (17) Av (18) sigma
%---------------------------------------------------------------------------------------%

    clear; clc; close all

%% Configurations
% soc range 
    soc_vec = 30; %0:10:100;
    % fitting할 soc 범위 지정, 0~100%, 10% 단위
    % 단일 SOC에 대해 진행할 경우 단일 SOC값만 입력. 
    
    save_path = 'G:\공유 드라이브\Battery Software Lab\Projects\LGES 2023\발표 및 공유자료\최종 보고서\Half_L\Natural'; %파라미터 결과 저장 폴더 경로 지정 
    save_check = 1; % 피규어 및 파라미터 데이터 저장 유무 선택, 0이면 결과 저장 x, 1이면 결과 저장 O
    % 주의: 동일한 폴더에 동일한 type_acf, type_dist 사용시 기존 파일 삭제 후 저장됨

% Fitting configuration
    type_weight = 1; % 0 for absolute error, 1 for relative error
    type_acf = 1; % 1 for anode, 2 for cathode, 3 for full cell (현재 구현되지 않는 상태)
    type_dist = 2; % 0 for DRT, 1 for DDT, 2 for integrated
    type_anode = 1; % 0 for base, 1 for natural 2 for blend
    num_iter = 100; %P2D 최적화 과정 최대 반복 횟수
    num_iter_dist = 10; % Dist 최적화 과정 최대 반복 횟수

% EIS data path
    if type_anode == 0
        % for base case
        path_folder = 'G:\공유 드라이브\BSL-Data\LGES\LG raw data\12_6cm2_soc10_EIS # Sample 1';%불러올 데이터 폴더 경로 지정
        path_file = 'PEIS_C09_%s_cycle_soc%d.csv'; % %s: anode or cathode, %d: soc percent
    
        % path_folder = 'G:\공유 드라이브\BSL-Data\LGES\전해질 DOE\전해액 변경_DOE 공유\processed_data';%불러올 데이터 폴더 경로 지정
        % path_file = '%s_천연흑연100_전해액변경 set_cell no_3_soc%d.csv'; % %s: anode or cathode, %d: soc percent
    
    elseif type_anode == 1
        % for anode Natural
        path_folder = 'G:\공유 드라이브\BSL-Data\LGES\Modeling data_2차 음극 2종\processed_data2';%불러올 데이터 폴더 경로 지정
        path_file = '%s_천역흑연100_음극_#2_soc%d.csv'; % %s: anode or cathode, %d: soc percent
    
    elseif type_anode == 2
        % for anode blending
        path_folder = 'G:\공유 드라이브\BSL-Data\LGES\Modeling data_2차 음극 2종\processed_data2';%불러올 데이터 폴더 경로 지정
        path_file = '%s_인조_천연 Blending음극_#1_soc%d.csv'; % %s: anode or cathode, %d: soc percent
    
    
    end 
%-----------------------------이 아래로는 수정 불필요-------------------------%

    if type_dist == 0
        dist = 'DRT'
    elseif type_dist == 1 
        dist = 'DDT'
    elseif type_dist == 2
        dist = 'DRT + DDT'
    end
    disp(dist);

    if type_acf == 1
        cell_type = 'Anode';
    elseif type_acf == 2
        cell_type = 'Cathode';
    elseif type_acf == 3
        cell_type = 'Full'; %현재 코드에선 구현 X
    end 

% Optimization options
    options= optimset('Display','iter','MaxIter',num_iter,'MaxFunEvals',1e5,...
        'TolFun',1e-10,'TolX',1e-10,'FinDiffType','central');
    options_dist= optimset('Display','iter','MaxIter',num_iter_dist,'MaxFunEvals',1e5,...
        'TolFun',1e-10,'TolX',1e-10,'FinDiffType','central');
% Parameters 
    bounds = [...
         0.05 20 % (1) R_itsc
         0.1 50; % (2) i0
         0.1 50; % (3) C_dl
         0.1 50; % (4) Ds
         0.1 50; % (5) kappa_el
         0.1 50; % (6) D_el
         0.1 50; % (7) Av
         0.1 10; % (8) Inductance (Henry)
         0.1 10; % (9) y_shift
         ]; 
    lb = bounds(:,1)*0.00001;
    ub = bounds(:,2)*10000;

    factors_ini = [1 1 1 1 1 1 1 1 1];

%Fitting start
para_sum = [];
for i = 1:length(soc_vec)

    SOC = soc_vec(i);
    name_file = sprintf(path_file, cell_type, SOC); 

% SOC and T (for initial guess - they are functions of soc, T)
    soc = SOC*0.01; % [1]
    T = 298.15; %[K]

%% Load and Pre-processing Data

    % load EIS data
    data = readmatrix([path_folder filesep name_file]);
    f_data = data(:,1);
    z_re_data = data(:,2);
    z_im_data = data(:,3);

    if z_im_data(1) < 0 
    z_im_data = - z_im_data;
    end

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
    % f_data = f_data(z_im_data<=0);
    % z_re_data = z_re_data(z_im_data<=0);
    % z_im_data = z_im_data(z_im_data<=0);
    % z_data = [z_re_data z_im_data];
    figure(1)
    plot(z_re_data,-z_im_data,'o'); hold on; grid on

    hold off;
%% Define weighting vector
    
    if type_weight == 1 % if minimizing the relative error
    weight = (z_re_data.^2 + z_im_data.^2).^(-0.5);
    weight_matrix = [weight weight];    
    elseif type_weight == 0 % if minimizing the absolute error
    weight_matrix = ones(size(z_re_data));
    end



fprintf('Start P2D fitting \n')

%% Call EIS model
   weighted_model = @(factors,f_data)BSL_func_EISmodel_V1_half_L(f_data, factors,soc,T,type_acf)...
       .*weight_matrix;
   weighted_data = z_data.*weight_matrix;

   tic;
   [factors_hat,resnorm,residual,~,~,~,jacobian_hat] ...
        = lsqcurvefit(weighted_model,factors_ini,...
                      f_data,weighted_data, lb, ub, options);
   toc;


%% Plot Results
[z_model0, paras0] = BSL_func_EISmodel_V1_half_L(f_data,factors_ini,soc,T,type_acf);
[z_model1, paras1] = BSL_func_EISmodel_V1_half_L(f_data,factors_hat,soc,T,type_acf);

fprintf('P2D fitting has successfully completed \n')
% Nyquist Plot
figure(2)
t = tiledlayout(1,2,"TileSpacing","compact",'Padding','compact');
nexttile
plot(z_data(:,1),-z_data(:,2),'ok','linewidth',1); hold on
%plot(z_model0(:,1),-z_model0(:,2),'ob','linewidth',1)
plot(z_model1(:,1),-z_model1(:,2),'or','linewidth',1)
%legend('Exp Data','Model Initial','Model Fit')
daspect ([1 1 2])

hold off;
    axis_limit = 1.1*max(max(abs(z_data)));
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
%    hold off
    grid on;    
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')
    legend('Exp Data','P2D')

% Zoom-in semicircle
nexttile
plot(z_data(:,1),-z_data(:,2),'ok','linewidth',1); hold on
%plot(z_model0(:,1),-z_model0(:,2),'ob','linewidth',1)
plot(z_model1(:,1),-z_model1(:,2),'or','linewidth',1)
%legend('Exp Data','Model Initial','Model Fit')


f_zoom_lb = 10; %[Hz] 
idx_zoom = f_data>f_zoom_lb;
axis_limit = 1.1*max(max(abs(z_data(idx_zoom,:))));
set(gca,'Box','on',... %Axis Properties: BOX   
'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
 hold off
 grid on;
 xlabel('Z_{re} [Ohm]')
 ylabel('-Z_{im} [Ohm]')
 legend('Exp Data','P2D')

set(gcf,'position',[400 600 1000 500])
pause(0.1) %for check P2D figure

%% Fitting Improvement by Distributed models.
fprintf(['Start P2D + ' dist ' fitting \n'])
% [dist] Initial guess and bounds
std_ini = 0.5;
factors_ini = [factors_hat std_ini std_ini]; % drt, ddt 분포 각각 고려

ub = factors_ini*10;
lb = factors_ini*0.1;

ub(end-1:end) = 5; %DRT&DDT Std
lb(end-1:end) = 0.01; %DRT&DDT Std

%%  [dist] Call EIS model
   weighted_model_dist = @(factors,f_data)BSL_func_EISmodel_V_half_Dist_integrated_L(f_data, factors,soc,T,type_acf,type_dist)...
       .*weight_matrix;
   weighted_data = z_data.*weight_matrix;

   tic;
   [factors_hat_dist, resnorm,residual,~,~,~,jacobian_hat] ...
        = lsqcurvefit(weighted_model_dist,factors_ini,...
                      f_data,weighted_data, lb, ub, options_dist);
   toc;

 %% [dist] Plot Results
[z_model2, paras2] = BSL_func_EISmodel_V_half_Dist_integrated_L(f_data,factors_hat_dist,soc,T,type_acf,type_dist);
fprintf(['P2D + ' dist ' fitting has successfully completed \n'])
% Nyquist Plot
figure(4)
t = tiledlayout(1,2,"TileSpacing","compact",'Padding','compact');
nexttile
plot(z_data(:,1),-z_data(:,2),'ok','linewidth',1); hold on
plot(z_model1(:,1),-z_model1(:,2),'or','linewidth',1)
plot(z_model2(:,1),-z_model2(:,2),'o','linewidth',1,'Color',[0.3010 0.7450 0.9330])

legend('Exp Data','P2D',['P2D' ' + ' dist])
axis_limit = 1.1*max(max(abs(z_data)));
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])

    grid on;
    title([cell_type ' soc' num2str(SOC) ' P2D' ' + ' dist])
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')
hold off


% Zoom-in semicircle
nexttile
plot(z_data(:,1),-z_data(:,2),'ok','linewidth',1); hold on
plot(z_model1(:,1),-z_model1(:,2),'or','linewidth',1)
plot(z_model2(:,1),-z_model2(:,2),'o','linewidth',1,'Color',[0.3010 0.7450 0.9330])
legend('Exp Data','P2D',['P2D' ' + ' dist])
 f_zoom_lb = 10; %[Hz] 
    idx_zoom = f_data>f_zoom_lb;
    axis_limit = 1.1*max(max(abs(z_data(idx_zoom,:))));
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
 
    grid on;
    title([cell_type ' soc' num2str(SOC) ' P2D' ' + ' dist ' Zoom-in'])
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')
hold off;

set(gcf,'position',[1400 600 1000 500])
end
