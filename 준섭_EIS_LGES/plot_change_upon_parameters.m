% This code is for comparison the shape changes upon each parameters

clc, clear, close all;


%% Configurations
% soc range 
    soc_vec = 50;

% EIS data path
    path_folder = 'G:\공유 드라이브\BSL-Data\LGES\LG raw data\12_6cm2_soc10_EIS # Sample 1';%불러올 데이터 폴더 경로 지정
    path_file = 'PEIS_C09_%s_cycle_soc%d.csv'; % %s 및 soc%d 부분은 고정, 그 외의 파일 이름만 수동으로 지정 ex) 파일 이름: soc20_cathode_LGES 라면, soc%d_cathode_LGES로 지정
    %주의: 데이터 형식 변경되면 적용 불가, 파일 이름에 % 특수기호 삭제 후 진행

% Fitting configuration
    type_weight = 1; % 0 for absolute error, 1 for relative error
    type_acf = 1; % 1 for anode, 2 for cathode
    soc_integ = 0;
    num_iter = 100; %P2D 최적화 과정 최대 반복 횟수
% Optimization options
    options= optimset('Display','iter','MaxIter',num_iter,'MaxFunEvals',1e5,...
        'TolFun',1e-8,'TolX',1e-8,'FinDiffType','central');


 if type_acf == 1
        cell_type = 'Anode';
    elseif type_acf == 2
        cell_type = 'Cathode';
    elseif type_acf == 3
        cell_type = 'Full'; %현재 코드에선 구현 X
 end 

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

    SOC = soc_vec;
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

    %% Define weighting vector
    
    if type_weight == 1 % if minimizing the relative error
    weight = (z_re_data.^2 + z_im_data.^2).^(-0.5);
    weight_matrix = [weight weight];    
    elseif type_weight == 0 % if minimizing the absolute error
    weight_matrix = ones(size(z_re_data));
    end



fprintf('Start P2D fitting \n')

%% Call EIS model
   weighted_model = @(factors,f_data)BSL_func_EISmodel_V1_half(f_data, factors,soc,T,type_acf,soc_integ)...
       .*weight_matrix;
   weighted_data = z_data.*weight_matrix;

   tic;
   [factors_hat,resnorm,residual,~,~,~,jacobian_hat] ...
        = lsqcurvefit(weighted_model,factors_ini,...
                      f_data,weighted_data, lb, ub, options);
   toc;


%% plot upon parameters
   para_name = {'i0' 'C_dl' 'Ds' 'K_el' 'D_el' 'Av'};

   figure
   t = tiledlayout(2,3,"TileSpacing","compact","Padding","compact");

   for i = 2:length(factors_hat)
       nexttile
       hold on;

       factor_vars = {};
       for k = 1:9
       factors = factors_hat;
       
       if i == 3
          factors(i) = factors_hat(i)*((4*10) - (4*k));
       else
          factors(i) = factors_hat(i)*0.2*k;
       end 
       
       factor_cell = {num2str(factors(i))};
       factor_vars(k) = factor_cell;

       [z_model, paras1] = BSL_func_EISmodel_V1_half(f_data,factors,soc,T,type_acf,0);
       
       plot(z_model(:,1),-z_model(:,2),'Marker','o')
       % pause(0.1);
       end 
       hold off;
       if i == 2 || i == 4 || i == 5;
            axis_limit = 1.1*max(max(abs(z_model)));
       elseif i == 3
           axis_limit = 1/10*max(max(abs(z_model)));
       elseif i == 6 
           axis_limit = 1.6*max(max(abs(z_model)));
       elseif i == 7
           axis_limit = 3.4*max(max(abs(z_model)));
       end
        set(gca,'Box','on',... %Axis Properties: BOX   
        'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
        'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
        'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    %    hold off
        grid on;    
        xlabel('Z_{re} [Ohm]')
        ylabel('-Z_{im} [Ohm]')
       title(para_name(i-1))
       legend(factor_vars, "Location","best")
       set(gcf,"Position",[100 100 1200 800])

   end
