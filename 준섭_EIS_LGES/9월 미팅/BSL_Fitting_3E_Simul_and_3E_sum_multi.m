%% 
% This function 
% (1) finds the best-fit parameters for a halfcells
% calls a EIS model function

% Version
% V1_3E_combined: EIS_function_V1_3E_combined: cathode-anode combined fitting
% more outputs: func, data, plots

% this code stands for 3E_sum & 3E2full multi plot

clear; clc; close all

addpath 'C:\Users\admin\Documents\GitHub\BSL_LGES2024'
%% Configurations
path_folder = 'G:\공유 드라이브\BSL-Data\LGES\LG raw data\12_6cm2_soc10_EIS # Sample 1';
    save_path = 'C:\Users\admin\Documents\GitHub\JunSub\준섭_EIS_LGES\3Esum & 3E2full';

    para_sum3E_sum = [];
    para_sum3E2full = [];

soc_vec = 70;
for i = 1:length(soc_vec)

SOC = soc_vec(i);
% EIS data path
    %path_folder = 'G:\Shared drives\BSL-Data\LGES\LG raw data\12_6cm2_soc10_EIS # Sample 2';
    name_file_n = sprintf('PEIS_C09_anode_cycle_soc%d.csv',SOC); %양극, 음극 각각 지정
    name_file_p = sprintf('PEIS_C09_cathode_cycle_soc%d.csv',SOC);


% SOC and T (for initial guess - they are functions of soc, T)
    soc = SOC*0.01; % [1]
    T = 298.15; %[K]

% Fitting configuration
    type_weight = 1; % 0 for absolute error, 1 for relative error %가중 방식 지정
    type_acf =3; % 1 for anode, 2 for cathode, 3 for full cell %셀 타입 지정 (3E 는 3으로 지정)

% Optimization options
    options= optimset('display','iter','MaxIter',100,'MaxFunEvals',1e5,...
        'TolFun',1e-8,'TolX',1e-8,'FinDiffType','central'); %옵션 설정, 반복 나타내고, 최대 반복횟수 100, step 및 function tolerance 설정. 중앙차분방식 미분. 최신 최적화 함수에는 optimoptions 사용. 

% Parameters 
    %   1         2     3     4      5     6         7      8     9   10    11     12 
    % R_itsc_n, i0_n, Cdl_n, Ds_n, Av_n, R_itsc_p, i0_p, Cdl_p, Ds_p, Av_p, k_el, D_el
    
    bounds = [...
         0.5 2 % (1) R_n
         0.1 500; % (2) i0_n
         0.1 500; % (3) C_dl_n
         0.1 50; % (4) Ds_n
         0.02 10; % (5) Av_n
         0.5 2;  % (6) R_p
         0.1 50; % (7) i0_p
         0.1 50; % (8) C_dl_p
         0.1 50; % (9) Ds_p
         0.1 10; % (10) Av_p
         0.1 10; % (11) kappa_el
         0.01 10; % (12) D_el
         ]; 
    lb = bounds(:,1); %바운더리 설정
    ub = bounds(:,2);

    factors_ini = ones(size(bounds,1),1); % 초기 파라미터를 factor라고 지정한 셈.
    

%% Structure Explained

    % Data and model func
    % for n - frequency points
    % w_data (n x 1) vector
    % [ Z_n_re Z_n_im Z_p_re Z_p_im] (n x 4) matrix

%% Load and Pre-processing Data

    % load EIS data
    data_n = load([path_folder filesep name_file_n]); %데이터 로드 (양,음극 각각)
    data_p = load([path_folder filesep name_file_p]);
    f_data = [data_n(:,1) data_p(:,1)]; %프리퀀시 데이터 통합
        if any(f_data(:,1)~=f_data(:,1)) %아마 2여야 하는데 잘못 적힌 듯 함
            error('freq range of EIS_p and EIS_n are different')
        end
    z_data = [data_n(:,2) data_n(:,3) data_p(:,2) data_p(:,3)]; %Z_real, imag 데이터 통합.
           % [ Z_n_re     Z_n_im      Z_p_re      Z_p_im]    (n x 4) matrix

           figure(1)
    plot(data_n(:,2),-data_n(:,3),'ro'); hold on; grid on %음극 초도 플랏
    plot(data_p(:,2),-data_p(:,3),'bo'); %양극 초도 플랏
    
    axis_limit = 1.1*max(max(abs(z_data)));%축 제한 값 설정 (최대값 1.1배로 설정)
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')

    % trim the high-frequency inductance part (positive imag data)
    ind_keep = data_n(:,3)<=0 | data_p(:,3)<=0; %0 이상인 imag data 데이터 자름. (0미만인 부분만 index 지정) 
    f_data = f_data(ind_keep); %프리퀀시도 트림
    z_data = z_data(ind_keep,:); %z_data 트림
    figure(1)
    plot(z_data(:,1),-z_data(:,2),'ro','MarkerFaceColor','red'); hold on; grid on %트림한 데이터 플랏 음극, 양극
    plot(z_data(:,3),-z_data(:,4),'bo','MarkerFaceColor','blue'); 

    

%% Define weighting vector
    
    if type_weight == 1 % if minimizing the relative error %에러 타입 별로 가중치 방식 설정
    weight_n = (z_data(:,1).^2 + z_data(:,2).^2).^(-0.5); %RMS 역수
    weight_p = (z_data(:,3).^2 + z_data(:,4).^2).^(-0.5);
    weight_matrix = [weight_n weight_n weight_p weight_p];    
    elseif type_weight == 0 % if minimizing the absolute error
    weight_matrix = ones(size(z_data));
    end




%% FITTING
%   Call EIS model
   weighted_model = @(factors,f_data)BSL_func_EISmodel_V1_3E(f_data, factors,soc,T,type_acf).*weight_matrix; %모델 예측값에 가중치 곱하여 모델 설정
   %model = @(factors,f_data)BSL_func_EISmodel_V1_3E(f_data, factors,soc,T,type_acf);
   weighted_data = z_data.*weight_matrix; %가중치 곱해진 데이터 설정
%   fitting
   tic;
   [factors_hat, resnorm,residual,~,~,~,jacobian_hat] ...
        = lsqcurvefit(weighted_model,factors_ini,...
                      f_data,weighted_data, lb, ub, options); %factor_ini 최적화. weighted model 및 data 사용. --> 초기점 설정 후 EIS_model 에서 계산해 파라미터 얻음
   toc;


%% Plot Results
[z_model0, paras0] = BSL_func_EISmodel_V1_3E(f_data,factors_ini,soc,T,type_acf);
[z_model1, paras1] = BSL_func_EISmodel_V1_3E(f_data,factors_hat,soc,T,type_acf);

z_data_sum0 = [z_model0(:,1)+z_model0(:,3) z_model0(:,2)+z_model0(:,4)];
z_data_sum1 = [z_model1(:,1)+z_model1(:,3) z_model1(:,2)+z_model1(:,4)];

% Nyquist Plot
    % line colors
    cmat_jet = jet(16); %컬러맵 배열 생성
    cmat_n = [cmat_jet(16,:);cmat_jet(14,:);cmat_jet(11,:)]; %애노드, 캐소드 각각 색 지정
    cmat_p = [cmat_jet(1,:);cmat_jet(3,:);cmat_jet(6,:)];

figure(2)
t = tiledlayout(1,2); %2x2 피규어 생성, set handle at t
nexttile
hold on
plot(z_data(:,1)+z_data(:,3),-[z_data(:,2)+z_data(:,4)],'o','linewidth',1,'color',cmat_n(1,:));  
% plot(z_data_sum0(:,1),-z_data_sum0(:,2),'o','linewidth',1,'color',cmat_n(3,:))
plot(z_data_sum1(:,1),-z_data_sum1(:,2),'o','linewidth',1,'color',cmat_n(2,:)) 
legend('Exp Data','3E sum')
%legend('Exp Data','Model Fit')
daspect ([1 1 2])


axis_limit = 1.1*max(max(abs(z_data)));
set(gca,'Box','on',... %Axis Properties: BOX   
'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
grid on;
hold off
title(sprintf('3Esum plot at soc%d',SOC))
xlabel('Z_{re} [Ohm]')
ylabel('-Z_{im} [Ohm]')

%Zoom-in semicircle
nexttile
hold on
plot(z_data(:,1)+z_data(:,3),-[z_data(:,2)+z_data(:,4)],'o','linewidth',1,'color',cmat_n(1,:));  
% plot(z_data_sum0(:,1),-z_data_sum0(:,2),'o','linewidth',1,'color',cmat_n(3,:))
plot(z_data_sum1(:,1),-z_data_sum1(:,2),'o','linewidth',1,'color',cmat_n(2,:)) 
legend('Exp Data','3E sum')
%legend('Exp Data','Model Fit')
daspect ([1 1 2])


f_zoom_lb = 10; %[Hz] 
idx_zoom = f_data>f_zoom_lb;
axis_limit = 1.1*max(max(abs(z_data(idx_zoom,:))));
set(gca,'Box','on',... %Axis Properties: BOX   
'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
grid on;
hold off
title(sprintf('3Esum plot at soc%d',SOC))
xlabel('Z_{re} [Ohm]')
ylabel('-Z_{im} [Ohm]')

t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'position',[100 100 1200 600])

para_sum3E_sum(1,1+SOC/10) = SOC; %파라미터 하나로 뭉치기
para_sum3E_sum(2:13,1+SOC/10) = paras1;
para_sum3E_sum(15,1+SOC/10) = resnorm;

fig_name = sprintf('3Esum plot at soc%d',SOC); %피규어 이름
% savefig(fullfile(save_path,fig_name));

%-------------------------------------------------------3E2full start

% Fitting configuration
    type_weight = 1; % 0 for absolute error, 1 for relative error %가중 방식 지정
    type_acf =4; % 1 for anode, 2 for cathode, 3 for 3E, 4 for full cell %셀 타입 지정 (3E 는 3으로 지정)


    

%% Structure Explained

    % Data and model func
    % for n - frequency points
    % w_data (n x 1) vector
    % [ Z_n_re Z_n_im Z_p_re Z_p_im] (n x 4) matrix

%% Load and Pre-processing Data

    % load EIS data
    data_n = load([path_folder filesep name_file_n]); %데이터 로드 (양,음극 각각)
    data_p = load([path_folder filesep name_file_p]);
    f_data = [data_n(:,1) data_p(:,1)]; %프리퀀시 데이터 통합
        if any(f_data(:,1)~=f_data(:,2)) %아마 2여야 하는데 잘못 적힌 듯 함
            error('freq range of EIS_p and EIS_n are different')
        end
    z_data = [data_n(:,2) data_n(:,3) data_p(:,2) data_p(:,3)]; %Z_real, imag 데이터 통합.
           % [ Z_n_re     Z_n_im      Z_p_re      Z_p_im]    (n x 4) matrix

        if type_acf == 4
           z_data = [z_data(:,1)+z_data(:,3) z_data(:,2)+z_data(:,4)]; % full cell의 경우 데이터 전부 통합
        end
        
    figure(3)
    plot(z_data(:,1),-z_data(:,2),'ro'); hold on; grid on %full cell 초도 플랏 (separator 값 두 번 들어갈 것으로 예상)
    
    axis_limit = 1.1*max(max(abs(z_data)));%축 제한 값 설정 (최대값 1.1배로 설정)
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')

    % trim the high-frequency inductance part (positive imag data)
    ind_keep = data_n(:,3)<=0 | data_p(:,3)<=0; %0 이상인 imag data 데이터 자름. (0미만인 부분만 index 지정) 
    f_data = f_data(ind_keep); %프리퀀시도 트림
    z_data = z_data(ind_keep,:); %z_data 트림
  
    

    

%% Define weighting vector
    
    if type_weight == 1 % if minimizing the relative error %에러 타입 별로 가중치 방식 설정
    weight_f = (z_data(:,1).^2 + z_data(:,2).^2).^(-0.5); %RMS 역수
    weight_matrix = [weight_f weight_f];    
    elseif type_weight == 0 % if minimizing the absolute error
    weight_matrix = ones(size(z_data));
    end




%% FITTING
%   Call EIS model
   weighted_model = @(factors,f_data)BSL_func_EISmodel_V1_3E_Sum(f_data, factors,soc,T,type_acf).*weight_matrix; %모델 예측값에 가중치 곱하여 모델 설정
   %model = @(factors,f_data)BSL_func_EISmodel_V1_3E(f_data, factors,soc,T,type_acf);
   weighted_data = z_data.*weight_matrix; %가중치 곱해진 데이터 설정
%   fitting
   tic;
   [factors_hat, resnorm,residual,~,~,~,jacobian_hat] ...
        = lsqcurvefit(weighted_model,factors_ini,...
                      f_data,weighted_data, lb, ub, options); %factor_ini 최적화. weighted model 및 data 사용. --> 초기점 설정 후 EIS_model 에서 계산해 파라미터 얻음
   toc;


%% Plot Results
[z_model2, paras2] = BSL_func_EISmodel_V1_3E_Sum(f_data,factors_hat,soc,T,type_acf);

% Nyquist Plot
    % line colors
    cmat_jet = jet(16); %컬러맵 배열 생성
    cmat_n = [cmat_jet(16,:);cmat_jet(14,:);cmat_jet(11,:)]; %애노드, 캐소드 각각 색 지정
    cmat_p = [cmat_jet(1,:);cmat_jet(3,:);cmat_jet(6,:)];

figure(4)
t = tiledlayout(1,2); %2x2 피규어 생성, set handle at t
nexttile %행따라 진행
plot(z_data(:,1),-z_data(:,2),'o','linewidth',1,'color',cmat_n(1,:)); hold on %full exp plot
% plot(z_model0(:,1),-z_model0(:,2),'o','linewidth',1,'color',cmat_n(3,:)) %full model ini plot
plot(z_data_sum1(:,1),-z_data_sum1(:,2),'bo','linewidth',1) % 3Esum
plot(z_model2(:,1),-z_model2(:,2),'o','linewidth',1,'color',cmat_n(2,:)) %3E2full model fit plot
legend('Exp Data','3E sum','3E2full')
daspect ([1 1 2])


    axis_limit = 1.1*max(max(abs(z_data)));
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    grid on;
    hold off
    title(sprintf('3E2full & 3Esum plot at soc%d',SOC))
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')

% Zoom-in semicircle
nexttile
plot(z_data(:,1),-z_data(:,2),'o','linewidth',1,'color',cmat_n(1,:)); hold on
% plot(z_model0(:,1),-z_model0(:,2),'o','linewidth',1,'color',cmat_n(3,:))
plot(z_data_sum1(:,1),-z_data_sum1(:,2),'bo','linewidth',1) % 3Esum
plot(z_model2(:,1),-z_model2(:,2),'o','linewidth',1,'color',cmat_n(2,:))
legend('Exp Data','3E sum','3E2full')



    f_zoom_lb = 10; %[Hz] 
    idx_zoom = f_data>f_zoom_lb;
    axis_limit = 1.1*max(max(abs(z_data(idx_zoom,:))));
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    grid on;
    hold off
    title(sprintf('3E2full & 3Esum plot at soc%d',SOC))
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')


t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'position',[100 100 1200 600])

fig_name = sprintf('3E2full & 3Esum plot at soc%d',SOC); %피규어 이름
% savefig(fullfile(save_path,fig_name));
%% Result Summary
   
Result.factors_hat = factors_hat;
Result.paras_hat = paras1;
Result.z_model = z_model1;

para_sum3E2full(1,1+SOC/10) = SOC; %파라미터 하나로 뭉치기
para_sum3E2full(2:13,1+SOC/10) = paras2;
para_sum3E2full(15,1+SOC/10) = resnorm;


% input('press enter to continue') 끊으며 확인하고 싶을 때 사용.
end

% save([save_path filesep  'para_sum3E_Simul.mat']','para_sum3E_sum');
% save([save_path filesep  'para_sum_Sum']','para_sum3E2full');

