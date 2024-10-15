%this code is for plotting distribution graph DRT or DDT ot both
clear; clc; close all

%% Configurations
% soc range
    soc_vec = 30; %0:10:100;
    % fitting할 soc 범위 지정, 0~100%, 10% 단위
    % 단일 SOC에 대해 진행할 경우 단일 SOC값만 입력. 
    % 이 코드는 단일 SOC에 대해서만 실행 or save 코드 추가
    
    save_path = 'C:\Users\admin\Documents\GitHub\JunSub\준섭_EIS_LGES\LGES_all_users'; %파라미터 결과 저장 폴더 경로 지정 
    save_check = 0; % 피규어 및 파라미터 데이터 저장 유무 선택, 0이면 결과 저장 x, 1이면 결과 저장 O
    % 주의: 동일한 폴더에 동일한 type_acf, type_dist 사용시 기존 파일 삭제 후 저장됨

% EIS data path
    path_folder = 'G:\공유 드라이브\BSL-Data\LGES\LG raw data\12_6cm2_soc10_EIS # Sample 1';%불러올 데이터 폴더 경로 지정
    path_file = 'PEIS_C09_%s_cycle_soc%d.csv'; % %s 및 soc%d 부분은 고정, 그 외의 파일 이름만 수동으로 지정 ex) 파일 이름: soc20_cathode_LGES 라면, soc%d_cathode_LGES로 지정
    %주의: 데이터 형식 변경되면 적용 불가, 파일 이름에 % 특수기호 삭제 후 진행

% Fitting configuration
    type_weight = 1; % 0 for absolute error, 1 for relative error
    type_acf = 2; % 1 for anode, 2 for cathode, 3 for full cell (현재 구현되지 않는 상태)
    type_dist = 2; % 0 for DRT, 1 for DDT, 2 for integrated
    num_iter = 5; %최적화 과정 최대 반복 횟수

% Distribution visuallize setting
    dist_target = 0; %0 for DRT, 1 for DDT, 2 for both
    vec_range_min = 0.1;
    vec_range_max = 2.1;
    vec_interval = 0.3;
    fixed_sigma = 0.5; %set stationary sigma



    std_vec = vec_range_min:vec_interval:vec_range_max; %range of sigma - equally applied
    
%-----------------------------이 아래로는 수정 불필요-------------------------%

    if type_dist == 0
        dist = 'DRT'
    elseif type_dist == 1 
        dist = 'DDT'
    elseif type_dist == 2
        dist = 'DRT + DDT'
    end


    if type_acf == 1
        cell_type = 'Anode';
    elseif type_acf == 2
        cell_type = 'Cathode';
    elseif type_acf == 3
        cell_type = 'Full'; %현재 코드에선 구현 X
    end 

% Optimization options
    options= optimset('Display','iter','MaxIter',num_iter,'MaxFunEvals',1e5,...
        'TolFun',1e-8,'TolX',1e-8,'FinDiffType','central');


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

    factors_ini = [1 1 1 1 1 1 1];

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
    f_data = f_data(z_im_data<=0);
    z_re_data = z_re_data(z_im_data<=0);
    z_im_data = z_im_data(z_im_data<=0);
    z_data = [z_re_data z_im_data];
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
   weighted_model = @(factors,f_data)BSL_func_EISmodel_V1_half(f_data, factors,soc,T,type_acf)...
       .*weight_matrix;
   weighted_data = z_data.*weight_matrix;

   tic;
   [factors_hat,resnorm,residual,~,~,~,jacobian_hat] ...
        = lsqcurvefit(weighted_model,factors_ini,...
                      f_data,weighted_data, lb, ub, options);
   toc;

%% P2D + DRT + DDT fitting 
std_ini = 0.5;
    factors_ini = [factors_hat std_ini std_ini]; % drt, ddt 분포 각각 고려

    ub = factors_ini*10;
    lb = factors_ini*0.1;

%%  [dist] Call EIS model
   weighted_model_dist = @(factors,f_data)BSL_func_EISmodel_V_half_Dist_integrated(f_data, factors,soc,T,type_acf,type_dist)...
       .*weight_matrix;
   weighted_data = z_data.*weight_matrix;

   tic;
   [factors_hat_dist, resnorm,residual,~,~,~,jacobian_hat] ...
        = lsqcurvefit(weighted_model_dist,factors_ini,...
                      f_data,weighted_data, lb, ub, options);
   toc;

%% P2D + DRT + DDT demonstration upon Std
factors_dist = factors_hat_dist; % drt, ddt 분포 각각 고려
        
figure(2)
hold on;
 
for k = 1:length(std_vec)
    
         if dist_target == 0 
             std_ini_drt = std_vec(k); 
             std_ini_ddt = fixed_sigma; 
         elseif dist_target == 1
             std_ini_drt = fixed_sigma; 
             std_ini_ddt = std_vec(k);
         elseif dist_target == 2
             std_ini_drt = std_vec(k); 
             std_ini_ddt = std_vec(k); 
         end 

        factors_dist(8) = std_ini_drt;
        factors_dist(9) = std_ini_ddt;

        [z_model2, paras2] = BSL_func_EISmodel_V_half_Dist_integrated(f_data,factors_dist,soc,T,type_acf,type_dist);
        
        imp{k} = z_model2;
        if k == 1
        plot(z_data(:,1),-z_data(:,2),'ok','linewidth',1); 
        end 
        plot(z_model2(:,1),-z_model2(:,2),'-','linewidth',1)
        
        legend_txt{k+1} = sprintf('\\sigma_{drt} = %.2f, \\sigma_{ddt} = %.2f', std_ini_drt, std_ini_ddt);
        
      
            
end 
hold off;
  axis_limit = 1.1*max(max(abs(z_data)));
            set(gca,'Box','on',... %Axis Properties: BOX   
            'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
            'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
            'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
        
            grid on;
            title([cell_type ' soc' num2str(SOC) ' P2D' ' + ' dist])
            xlabel('Z_{re} [Ohm]')
            ylabel('-Z_{im} [Ohm]')
legend_txt{1} = 'Exp data';
legend('Exp Data', legend_txt)
set(gcf,'position',[300 300 1000 1000])

figure(3)
    hold on;
 
for k = 1:length(std_vec)
    % Zoom-in semicircle
         if dist_target == 0 
             std_ini_drt = std_vec(k); 
             std_ini_ddt = fixed_sigma; 
         elseif dist_target == 1
             std_ini_drt = fixed_sigma; 
             std_ini_ddt = std_vec(k);
         elseif dist_target == 2
             std_ini_drt = std_vec(k); 
             std_ini_ddt = std_vec(k); 
         end 

        factors_dist(8) = std_ini_drt;
        factors_dist(9) = std_ini_ddt;
        [z_model2, paras2] = BSL_func_EISmodel_V_half_Dist_integrated(f_data,factors_dist,soc,T,type_acf,type_dist);
        
        imp{k} = z_model2;
        if k == 1
        plot(z_data(:,1),-z_data(:,2),'ok','linewidth',1); 
        end 
        plot(z_model2(:,1),-z_model2(:,2),'-','linewidth',1)
        
  
end
 f_zoom_lb = 10; %[Hz] 
        idx_zoom = f_data>f_zoom_lb;
        axis_limit = 1.1*max(max(abs(z_data(idx_zoom,:))));
        set(gca,'Box','on',... %Axis Properties: BOX   
        'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
        'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
        'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
     
        grid on;
        title([cell_type ' soc' num2str(SOC) ' P2D' ' + ' dist  ' Zoom-in'])
        xlabel('Z_{re} [Ohm]')
        ylabel('-Z_{im} [Ohm]')
hold off;
legend('Exp Data', legend_txt)
set(gcf,'position',[1400 300 1000 1000])
end 