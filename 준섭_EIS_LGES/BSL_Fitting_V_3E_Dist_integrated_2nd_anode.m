%% 
% This function 
% (1) finds the best-fit parameters for a halfcells calls a EIS model function
% (2) improve fitting result by addapting additional distribution model

%Configuration 설정 후 최적화 진행
% Cylindrical data는 0:20:100 간격
% 3E data는 0:10:100 간격
%---------------------------------------------------------------------------------------%

    clear; clc; close all
% for n = 1:2 %for multilple acquisition of data

%% Configurations
    save_path = 'C:\Users\admin\Documents\GitHub\JunSub\준섭_EIS_LGES\10월 미팅\multi-soc_DRT_plot_2nd_anode_Natural'; %파라미터 결과 저장 폴더 경로 지정 
    save_check = 1; % 피규어 및 파라미터 데이터 저장 유무 선택, 0이면 결과 저장 x, 1이면 결과 저장 O
    % 주의: 동일한 폴더에 동일한 type_acf, type_dist 사용시 기존 파일 삭제 후 저장됨

% soc range 
    soc_start = 30;
    soc_end = 70;   
    soc_interval = 20;

    multi_start = 20;
    multi_end = 80;
% Fitting configuration
    type_weight = 1; % 0 for absolute error, 1 for relative error
    type_cell = 1; % 1 for 3E cell, 2 for cylindrical cell (only for 3E_Sum)
    type_acf = 4; % 0 for full, 1 for anode, 2 for cathode, 3 for 3E_Sum, 4 for 3E_Simul
    type_anode = 2; %1 for blending, 2 for natural graphite 100%
    type_dist = 0; % 0 for DRT, 1 for DDT, 2 for integrated
    soc_integ = 1; % 0 for inactive, 1 for active

    num_iter = 100; %P2D 최적화 과정 최대 반복 횟수
    num_iter_dist = 100; % Dist 최적화 과정 최대 반복 횟수

% Temperature
    T = 298.15; %[K]

% EIS data path
    if type_cell == 1    
        path_folder = 'G:\공유 드라이브\BSL-Data\LGES\Modeling data_2차 음극 2종';
        path_file = '%s_soc%d.txt';
    elseif type_cell == 2
        path_folder = 'G:\공유 드라이브\BSL-Data\LGES\2170_Cell_for modeling';
        path_file = '%s_soc%d.csv';
    end
    
    if type_cell == 2 && type_acf ~= 3
       error('Cylindrical full cell only can use 3E_Sum (full)')
    end

%-----------------------------이 아래로는 수정 불필요-------------------------%
% soc setting
    if type_cell == 1    
        soc_vec = soc_start:soc_interval:soc_end;  
    elseif type_cell == 2
        soc_vec = soc_start:20:soc_end;  
    end 
        multi_soc_range = multi_start:10:multi_end; % for soc integrate
% Naming configuration
    if type_acf == 0
        cell_type = 'Full'
    elseif type_acf == 1
        cell_type = 'Anode';
    elseif type_acf == 2
        cell_type = 'Cathode';
    elseif type_acf == 3
        cell_type = '3E sum';
    elseif type_acf == 4
        cell_type = '3E simul';
    end 
    

    if type_anode == 1
        if type_cell == 1
            anode_type = '인조_천연 Blending음극_#1';
        elseif type_cell == 2
            anode_type = 'Set2_인조_천연_Blending_Cell No_2';
        end 
    elseif type_anode == 2
        if type_cell == 1
            anode_type = '천역흑연100_음극_#2';
        elseif type_cell == 2
            anode_type = 'Set1_천연100_Cell No_2';
        end 
    end
    

    if type_anode == 1
        if type_cell == 1 
            name_type = '3E Art&Nat blending';
        elseif type_cell == 2
            name_type = 'Cylindrical Art&Nat blending';
        end

    elseif type_anode == 2
         if type_cell == 1 
            name_type = '3E Natural 100%';
         elseif type_cell == 2
            name_type = 'Cylindrical Natural 100%';
         end
    end 

    if soc_integ == 0
        if type_dist == 0
            dist = 'DRT';
        elseif type_dist == 1 
            dist = 'DDT';
        elseif type_dist == 2
            dist = 'DRT + DDT';
        end
    elseif soc_integ == 1
        if type_dist == 0
            dist = 'soc integrated + DRT';
        elseif type_dist == 1 
            dist = 'soc integrated + DDT';
        elseif type_dist == 2
            dist = 'soc integrated + DRT + DDT';
        end
    end 
    disp([cell_type ' ' dist]);

    %% half_cell script
if type_acf == 0 || type_acf == 1 || type_acf == 2
% Optimization options
    options= optimset('Display','iter','MaxIter',num_iter,'MaxFunEvals',1e5,...
        'TolFun',1e-8,'TolX',1e-8,'FinDiffType','central');

    options_dist= optimset('Display','iter','MaxIter',num_iter_dist,'MaxFunEvals',1e5,...
        'TolFun',1e-8,'TolX',1e-8,'FinDiffType','central');

    
%Fitting start
para_sum =zeros(19,11);
for i = 1:length(soc_vec)
    % Parameters 
    bounds = [...
         0.5 2 % (1) R_itsc
         0.1 50; % (2) i0
         0.1 50; % (3) C_dl
         0.001 50; % (4) Ds
         0.001 100; % (5) kappa_el
         0.001 100; % (6) D_el
         0.1 10; % (7) Av
         ]; 
    lb = 0.001*bounds(:,1);
    ub = 30*bounds(:,2);
    
    factors_ini = [1 1 1 1 1 1 1];
    SOC = soc_vec(i);
   
    name_file = sprintf(path_file, anode_type, SOC); 
% SOC and T (for initial guess - they are functions of soc, T)
    soc = SOC*0.01; % [1]
    
% SOC integration
if soc_integ == 1 && i == 1
   [f_data, z_data_integ, z_model_integ, z_model_dist_integ, paras_integ ,paras_integ_dist] = BSL_func_3E_calc_el_multi_soc_2nd_anode(path_folder, path_file, multi_soc_range, cell_type, anode_type);
    for j = 1:length(multi_soc_range)
        SOC = multi_soc_range(j);
        figure(j)
        t = tiledlayout(1,2,"TileSpacing","compact",'Padding','compact');
        nexttile
        plot(z_data_integ(:,2*j-1),-z_data_integ(:,2*j),'ok','linewidth',1); hold on
        plot(z_model_integ(:,2*j-1),-z_model_integ(:,2*j),'or','linewidth',1)
        plot(z_model_dist_integ(:,2*j-1),-z_model_dist_integ(:,2*j),'o','linewidth',1,'Color',[0.3010 0.7450 0.9330])
        
        legend('Exp Data','P2D',['P2D' ' + ' dist])
        axis_limit = 1.1*max(max(abs(z_data_integ(:,2*j-1:2*j))));
            set(gca,'Box','on',... %Axis Properties: BOX   
            'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
            'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
            'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
        
            grid on;
            title([cell_type ' soc' num2str(SOC) ' P2D'         ' + ' dist])
            xlabel('Z_{re} [Ohm]')
            ylabel('-Z_{im} [Ohm]')
        hold off
        
        
        % Zoom-in semicircle
        nexttile
        plot(z_data_integ(:,2*j-1),-z_data_integ(:,2*j),'ok','linewidth',1); hold on
        plot(z_model_integ(:,2*j-1),-z_model_integ(:,2*j),'or','linewidth',1)
        plot(z_model_dist_integ(:,2*j-1),-z_model_dist_integ(:,2*j),'o','linewidth',1,'Color',[0.3010 0.7450 0.9330])
        legend('Exp Data','P2D',['P2D' ' + ' dist])
         f_zoom_lb = 10; %[Hz] 
            idx_zoom = f_data>f_zoom_lb;
            axis_limit = 1.1*max(max(abs(z_data_integ(idx_zoom,:))));
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
        
        if save_check == 1
            fig_name = sprintf([anode_type ' ' cell_type ' soc_%d' ' P2D' ' + ' dist ' .fig'],SOC); %피규어 저장
            savefig(fullfile(save_path,fig_name));
        end
    end 
        

        if save_check == 1
            paras_integ(2:end+1,:) = paras_integ;
            paras_integ(1,:) = str2double(split(num2str(multi_soc_range)));
            save([save_path filesep anode_type ' ' cell_type  ' P2D' ' + ' dist ' para_sum_integ '],'paras_integ')

            paras_integ_dist(2:end+1,:) = paras_integ_dist;
            paras_integ_dist(1,:) = str2double(split(num2str(multi_soc_range)));
            save([save_path filesep anode_type ' ' cell_type  ' P2D' ' + ' dist ' para_sum_integ_dist '],'paras_integ_dist')

            fprintf('Data saved successfully.\n');
        else
            fprintf('Data not saved.\n')
        end
elseif soc_integ == 1 && i > 1
       Useless(i) = i;
elseif soc_integ == 0
%% Load and Pre-processing Data

    % load EIS data
        data = readmatrix([path_folder filesep name_file]);
        f_data = data(:,1);

        if type_acf == 0
           z_re_data = data(:,2);
           z_im_data = -data(:,3);
        elseif type_acf == 1
           z_re_data = -data(:,4);
           z_im_data = data(:,5);
        elseif type_acf == 2
           z_re_data = data(:,6);
           z_im_data = -data(:,7);
        end 

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
   weighted_model = @(factors,f_data)BSL_func_EISmodel_V1_half_2nd_anode(f_data, factors,soc,T,type_acf,type_anode)...
       .*weight_matrix;
   weighted_data = z_data.*weight_matrix;

   tic;
   [factors_hat,resnorm,residual,~,~,~,jacobian_hat] ...
        = lsqcurvefit(weighted_model,factors_ini,...
                      f_data,weighted_data, lb, ub, options);
   toc;


%% Plot Results
% [z_model0, paras0] = BSL_func_EISmodel_V1_half_2nd_anode(f_data,factors_ini,soc,T,type_acf,type_anode);
[z_model1, paras1] = BSL_func_EISmodel_V1_half_2nd_anode(f_data,factors_hat,soc,T,type_acf,type_anode);

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

    ub = factors_ini*30;
    lb = factors_ini*0.01;

%%  [dist] Call EIS model
   weighted_model_dist = @(factors,f_data)BSL_func_EISmodel_V_half_Dist_integrated_2nd_anode(f_data, factors,soc,T,type_acf,type_dist,type_anode)...
       .*weight_matrix;
   weighted_data = z_data.*weight_matrix;

   tic;
   [factors_hat_dist, resnorm,residual,~,~,~,jacobian_hat] ...
        = lsqcurvefit(weighted_model_dist,factors_ini,...
                      f_data,weighted_data, lb, ub, options_dist);
   toc;

 %% [dist] Plot Results
[z_model2, paras2] = BSL_func_EISmodel_V_half_Dist_integrated_2nd_anode(f_data,factors_hat_dist,soc,T,type_acf,type_dist,type_anode);
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
    title([name_type ' ' cell_type ' soc' num2str(SOC)])
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
    title([name_type ' ' cell_type ' soc' num2str(SOC) ' Zoom-in'])
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')
hold off;

set(gcf,'position',[1400 600 1000 500])



% %Beta comparison
% figure()
% hold on;
% plot(real(betaa_sum(1,:)),imag(betaa_sum(1,:)),'bo-');
% plot(real(betaa_sum(2,:)),imag(betaa_sum(2,:)),'ro-');
% hold off;
% legend('beta integ','beta eq');
% 
% %Zloc comparison
% zloca_sum = (betaa_sum).^-1;
% zlocc_sum = (betac_sum).^-1;
% 
% figure()
% hold on;
% plot(real(zloca_sum(1,:)),-imag(zloca_sum(1,:)),'bo-');
% plot(real(zloca_sum(2,:)),-imag(zloca_sum(2,:)),'ro-');
% hold off;
% legend('zloc integ','zloc eq');
%% Result Summary
   
Result.factors_hat = factors_hat;
Result.paras_hat = paras1;
Result.z_model = z_model1;


Result.factors_hat_dist = factors_hat_dist;
Result.paras_hat_dist = paras2;
Result.z_model_dist = z_model2;

Result_table = [[paras1;0;0] paras2];

if save_check == 1
fig_name = sprintf([anode_type ' ' cell_type ' soc_%d' ' P2D' ' + ' dist ' .fig'],SOC); %피규어 저장
savefig(fullfile(save_path,fig_name));
end

para_sum(1,1+SOC/10) = SOC;
para_sum(2:8,1+SOC/10) = paras1;
para_sum(11:19,1+SOC/10) = paras2;
% input('press enter to continue');

fprintf(['Process ' num2str(i) ' of ' num2str(length(soc_vec)) '\n']);
  
end
end

     if save_check == 1
            save([save_path filesep anode_type ' ' cell_type  ' P2D' ' + ' dist ' para_sum '],'para_sum')
            fprintf('Data saved successfully.\n');
        else
            fprintf('Data not saved.\n')
     end
%% 3E_Sum script
elseif type_acf == 3
    % Optimization options
    if type_cell == 1
        options= optimset('Display','iter','MaxIter',num_iter,'MaxFunEvals',1e5,...
            'TolFun',1e-8,'TolX',1e-8,'FinDiffType','central');
    
        options_dist= optimset('Display','iter','MaxIter',num_iter_dist,'MaxFunEvals',1e5,...
            'TolFun',1e-8,'TolX',1e-8,'FinDiffType','central');
    elseif type_cell == 2
        options= optimset('Display','iter','MaxIter',num_iter,'MaxFunEvals',1e5,...
        'TolFun',1e-9,'TolX',1e-9,'FinDiffType','central');
    
        options_dist= optimset('Display','iter','MaxIter',num_iter_dist,'MaxFunEvals',1e5,...
        'TolFun',1e-9,'TolX',1e-9,'FinDiffType','central'); 
    end 
   
    
    para_sum = zeros(32,11);
    for i = 1:length(soc_vec)
    
        % Parameters 
    %   1         2     3     4      5     6         7      8     9   10    11     12 
    % R_itsc_n, i0_n, Cdl_n, Ds_n, Av_n, R_itsc_p, i0_p, Cdl_p, Ds_p, Av_p, k_el, D_el

    if type_cell == 1
          bounds = [...
         0.05 2 % (1) R_n
         0.1 50; % (2) i0_n
         0.1 50; % (3) C_dl_n
         0.1 50; % (4) Ds_n
         0.02 10; % (5) Av_n
         0.05 2;  % (6) R_p
         0.1 50; % (7) i0_p
         0.1 2; % (8) C_dl_p
         0.1 50; % (9) Ds_p
         0.1 10; % (10) Av_p
         0.001 100; % (11) kappa_el
         0.001 100; % (12) D_el
         ]; 

        lb = 0.01*bounds(:,1); %바운더리 설정
        ub = 30*bounds(:,2);
        factors_ini = ones(1,size(bounds,1));

    elseif type_cell == 2
          bounds = [...
         0.05 2 % (1) R_n
         0.1 50; % (2) i0_n
         0.001 50; % (3) C_dl_n
         0.1 50; % (4) Ds_n
         0.02 10; % (5) Av_n
         0.05 2;  % (6) R_p
         0.1 50; % (7) i0_p
         0.001 50; % (8) C_dl_p
         0.1 50; % (9) Ds_p
         0.1 50; % (10) Av_p
         0.001 50; % (11) kappa_el
         0.001 50; % (12) D_el
         ]; 
        
        lb = 0.1*bounds(:,1); %바운더리 설정
        ub = 10*bounds(:,2);
        factors_ini = 0.01*ones(size(bounds,1),1); %
    end 
        
        SOC = soc_vec(i);

        soc = SOC*0.01; % [1]

        % soc integration
        if soc_integ == 1 && i == 1
           [f_data, z_data_integ, z_model_integ, z_model_dist_integ, paras_integ ,paras_integ_dist] = BSL_func_3E_calc_el_multi_soc_2nd_anode(path_folder, path_file, multi_soc_range, cell_type, anode_type);
            for j = 1:length(multi_soc_range)
                SOC = multi_soc_range(j);
                figure(j)
                t = tiledlayout(1,2,"TileSpacing","compact",'Padding','compact');
                nexttile
                plot(z_data_integ(:,2*j-1),-z_data_integ(:,2*j),'ok','linewidth',1); hold on
                plot(z_model_integ(:,2*j-1),-z_model_integ(:,2*j),'or','linewidth',1)
                plot(z_model_dist_integ(:,2*j-1),-z_model_dist_integ(:,2*j),'o','linewidth',1,'Color',[0.3010 0.7450 0.9330])
                
                legend('Exp Data','P2D',['P2D' ' + ' dist])
                axis_limit = 1.1*max(max(abs(z_data_integ(:,2*j-1:2*j))));
                    set(gca,'Box','on',... %Axis Properties: BOX   
                    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
                    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
                    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
                
                    grid on;
                    title([cell_type ' soc' num2str(SOC) ' P2D'         ' + ' dist])
                    xlabel('Z_{re} [Ohm]')
                    ylabel('-Z_{im} [Ohm]')
                hold off
                
                
                % Zoom-in semicircle
                nexttile
                plot(z_data_integ(:,2*j-1),-z_data_integ(:,2*j),'ok','linewidth',1); hold on
                plot(z_model_integ(:,2*j-1),-z_model_integ(:,2*j),'or','linewidth',1)
                plot(z_model_dist_integ(:,2*j-1),-z_model_dist_integ(:,2*j),'o','linewidth',1,'Color',[0.3010 0.7450 0.9330])
                legend('Exp Data','P2D',['P2D' ' + ' dist])
                 f_zoom_lb = 10; %[Hz] 
                    idx_zoom = f_data>f_zoom_lb;
                    axis_limit = 1.1*max(max(abs(z_data_integ(idx_zoom,:))));
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
                
                if save_check == 1
                    fig_name = sprintf([anode_type ' ' cell_type ' soc_%d' ' P2D' ' + ' dist ' .fig'],SOC); %피규어 저장
                    savefig(fullfile(save_path,fig_name));
                end
            end 
                
        
                if save_check == 1
                    paras_integ(2:end+1,:) = paras_integ;
                    paras_integ(1,:) = str2double(split(num2str(multi_soc_range)));
                    save([save_path filesep anode_type ' ' cell_type  ' P2D' ' + ' dist ' para_sum_integ '],'paras_integ')
        
                    paras_integ_dist(2:end+1,:) = paras_integ_dist;
                    paras_integ_dist(1,:) = str2double(split(num2str(multi_soc_range)));
                    save([save_path filesep anode_type ' ' cell_type  ' P2D' ' + ' dist ' para_sum_integ_dist '],'paras_integ_dist')
        
                    fprintf('Data saved successfully.\n');
                else
                    fprintf('Data not saved.\n')
                end
        elseif soc_integ == 1 && i > 1
               Useless(i) = i;
        elseif soc_integ == 0
        %% Structure Explained
            % Data and model func
            % for n - frequency points
            % w_data (n x 1) vector
            % [ Z_n_re Z_n_im Z_p_re Z_p_im] (n x 4) matrix
        
            % load EIS data
            name_file = sprintf(path_file, anode_type, SOC); 
            
            data = readmatrix([path_folder filesep name_file]);
            if type_cell == 1
                data_n = [data(:,1) -data(:,4) -data(:,5)];
                data_p = [data(:,1) data(:,6) data(:,7)];
            
                f_data = [data_n(:,1) data_p(:,1)]; %프리퀀시 데이터 통합
               
                
                    if any(f_data(:,1)~=f_data(:,2))
                        error('freq range of EIS_p and EIS_n are different')
                    end
            
                z_data = [data_n(:,2) data_n(:,3) data_p(:,2) data_p(:,3)]; %Z_real, imag 데이터 통합.
                       % [Z_n_re     Z_n_im      Z_p_re      Z_p_im]    (n x 4) matrix
                z_data = [z_data(:,1)+z_data(:,3) -(z_data(:,2)+z_data(:,4))];
                           %[Z_Real_sum Z_imag_sum]; 
           
            elseif type_cell == 2
                f_data = data(:,1);
                z_data(:,1) = data(:,2); 
                z_data(:,2) = data(:,3);
            end
            figure(1)
            plot(z_data(:,1),-z_data(:,2),'ro'); hold on; grid on %full cell 초도 플랏 (separator 값 두 번 들어갈 것으로 예상)
            
            axis_limit = 1.1*max(max(abs(z_data)));%축 제한 값 설정 (최대값 1.1배로 설정)
            set(gca,'Box','on',... %Axis Properties: BOX   
            'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
            'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
            'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
            hold off;
            xlabel('Z_{re} [Ohm]')
            ylabel('-Z_{im} [Ohm]')
            
            pause(0.01);
            % trim the high-frequency inductance part (positive imag data)
            if type_cell == 1
                ind_keep = data_n(:,3) >= 0 & data_p(:,3) >= 0; % imag 0 이상인 imag data 데이터 자름. (0미만인 부분만 index 지정) 
                f_data = f_data(ind_keep); %프리퀀시도 트림
                z_data = z_data(ind_keep,:); %z_data 트림
            elseif type_cell == 2
                ind_keep = data(:,3) <= 0; %0 이상인 imag data 데이터 자름. (0미만인 부분만 index 지정) 
                f_data = f_data(ind_keep);
                z_data = z_data(ind_keep,:); %z_data 트림
            end 
        %% Define weighting vector
    
            if type_weight == 1 % if minimizing the relative error %에러 타입 별로 가중치 방식 설정
            weight_f = (z_data(:,1).^2 + z_data(:,2).^2).^(-0.5); %RMS 역수
            weight_matrix = [weight_f weight_f];    
            elseif type_weight == 0 % if minimizing the absolute error
            weight_matrix = ones(size(z_data));
            end

        %% FITTING
            %   Call EIS model
               weighted_model = @(factors,f_data)BSL_func_EISmodel_V1_3E_2nd_anode(f_data, factors,soc,T,type_acf,type_anode,type_cell).*weight_matrix; %모델 예측값에 가중치 곱하여 모델 설정
               weighted_data = z_data.*weight_matrix; %가중치 곱해진 데이터 설정
            %   fitting
               tic;
               [factors_hat, resnorm,residual,~,~,~,jacobian_hat] ...
                    = lsqcurvefit(weighted_model,factors_ini,...
                                  f_data,weighted_data, lb, ub, options); %factor_ini 최적화. weighted model 및 data 사용. --> 초기점 설정 후 EIS_model 에서 계산해 파라미터 얻음
               toc;
            
            
        %% Plot Results
            [z_model1, paras1] = BSL_func_EISmodel_V1_3E_2nd_anode(f_data,factors_hat,soc,T,type_acf,type_anode,type_cell);
            
            % Nyquist Plot
                % line colors
                cmat_jet = jet(16); %컬러맵 배열 생성
                cmat_n = [cmat_jet(16,:);cmat_jet(14,:);cmat_jet(11,:)]; %애노드, 캐소드 각각 색 지정
                cmat_p = [cmat_jet(1,:);cmat_jet(3,:);cmat_jet(6,:)];
            
            figure(2)
            t = tiledlayout(1,2); %2x2 피규어 생성, set handle at t
            nexttile %행따라 진행
            plot(z_data(:,1),-z_data(:,2),'o','linewidth',1,'color',cmat_n(1,:)); hold on %full exp plot
            % plot(z_model0(:,1),-z_model0(:,2),'o','linewidth',1,'color',cmat_n(3,:)) %full model ini plot
            plot(z_model1(:,1),-z_model1(:,2),'o','linewidth',1,'color',cmat_n(2,:)) %3E_cum model fit plot
            legend('Exp Data','3E sum')
            daspect ([1 1 2])

            axis_limit = 1.1*max(max(abs(z_data)));
            set(gca,'Box','on',... %Axis Properties: BOX   
            'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
            'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
            'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
            grid on;
            hold off
           
            title(sprintf('%s 3E sum plot at soc%d',name_type,SOC))
            xlabel('Z_{re} [Ohm]')
            ylabel('-Z_{im} [Ohm]')

            % Zoom-in semicircle
            nexttile
            plot(z_data(:,1),-z_data(:,2),'o','linewidth',1,'color',cmat_n(1,:)); hold on
            % plot(z_model0(:,1),-z_model0(:,2),'o','linewidth',1,'color',cmat_n(3,:))
            plot(z_model1(:,1),-z_model1(:,2),'o','linewidth',1,'color',cmat_n(2,:))
            legend('Exp Data','3E sum')

            f_zoom_lb = 10; %[Hz] 
            idx_zoom = f_data>f_zoom_lb;
            axis_limit = 1.1*max(max(abs(z_data(idx_zoom,:))));
            set(gca,'Box','on',... %Axis Properties: BOX   
            'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
            'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
            'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
            if type_cell == 2
                xlim([0.97*z_data(1,1) 1.02/1.1*axis_limit])
                ylim([0 1.02/1.1*axis_limit - 0.97*z_data(1,1)])
            end 
            grid on;
            hold off
            title(sprintf('%s 3E sum plot at soc%d',name_type,SOC))
            xlabel('Z_{re} [Ohm]')
            ylabel('-Z_{im} [Ohm]')
        
        
            t.TileSpacing = 'compact';
            t.Padding = 'compact';
            set(gcf,'position',[100 100 1200 600])

            pause(0.1); %for showing P2D figure
         %% Fitting Improvement by Distributed models.
            std_ini = 0.5;    
            factors_ini = [factors_hat std_ini std_ini std_ini std_ini]; % drta,c / ddt,ac 분포 각각 고려
           
            if type_cell == 1
                ub = factors_ini*100;
                lb = factors_ini*0.01;
            elseif type_cell == 2
                ub = [ub 5 5 5 5];
                lb = [lb 0.05 0.05 0.05 0.05];
            end 

         %% [dist] Call EIS model
            weighted_model_dist = @(factors,f_data)BSL_func_EISmodel_V_3E_Dist_integrated_2nd_anode(f_data, factors,soc,T,type_acf,type_dist,type_anode,type_cell)...
               .*weight_matrix;
            weighted_data = z_data.*weight_matrix;
        
            tic;
            [factors_hat_dist, resnorm,residual,~,~,~,jacobian_hat] ...
                = lsqcurvefit(weighted_model_dist,factors_ini,...
                              f_data,weighted_data, lb, ub, options_dist);
            toc;

            %% [dist] fitting result
            % Nyquist Plot
            [z_model2, paras2] = BSL_func_EISmodel_V_3E_Dist_integrated_2nd_anode(f_data,factors_hat_dist,soc,T,type_acf,type_dist,type_anode,type_cell);
            figure(3)
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
            title([name_type ' ' cell_type ' soc' num2str(SOC)])
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
         
            if type_cell == 2
                xlim([0.97*z_data(1,1) 1.02/1.1*axis_limit])
                ylim([0 1.02/1.1*axis_limit - 0.97*z_data(1,1)])
            end 
            grid on;
            title([name_type ' ' cell_type ' soc' num2str(SOC) ' Zoom-in'])
            xlabel('Z_{re} [Ohm]')
            ylabel('-Z_{im} [Ohm]')
            hold off;
            
            set(gcf,'position',[1400 600 1000 500])

        %% Result save
            para_sum(1,1+SOC/10) = SOC; %파라미터 하나로 뭉치기
            para_sum(2:13,1+SOC/10) = paras1;
            para_sum(15:30,1+SOC/10) = paras2;
            para_sum(32,1+SOC/10) = resnorm; 

            if save_check == 1
            fig_name = sprintf([anode_type ' ' cell_type ' soc_%d' ' P2D' ' + ' dist ' .fig'],SOC); %피규어 저장
            savefig(fullfile(save_path,fig_name));
            end
            
    end 
    end
            
            if save_check == 1
                save([save_path filesep anode_type ' ' cell_type  ' P2D' ' + ' dist ' para_sum '],'para_sum')
                fprintf('Data saved successfully.\n');
            else
                fprintf('Data not saved.\n')
            end
%% 3E_Simul script
elseif type_acf == 4 
    % Optimization options
    options= optimset('Display','iter','MaxIter',num_iter,'MaxFunEvals',1e5,...
        'TolFun',1e-8,'TolX',1e-8,'FinDiffType','central');

    options_dist= optimset('Display','iter','MaxIter',num_iter_dist,'MaxFunEvals',1e5,...
        'TolFun',1e-8,'TolX',1e-8,'FinDiffType','central');
    
    para_sum = zeros(32,11);
    for i = 1:length(soc_vec)

        % Parameters 
    %   1         2     3     4      5     6         7      8     9   10    11     12 
    % R_itsc_n, i0_n, Cdl_n, Ds_n, Av_n, R_itsc_p, i0_p, Cdl_p, Ds_p, Av_p, k_el, D_el
    
    bounds = [...
         0.5 2 % (1) R_n
         0.1 50; % (2) i0_n
         0.1 50; % (3) C_dl_n
         0.1 50; % (4) Ds_n
         0.02 10; % (5) Av_n
         0.5 2;  % (6) R_p
         0.1 50; % (7) i0_p
         0.1 50; % (8) C_dl_p
         0.1 50; % (9) Ds_p
         0.1 10; % (10) Av_p
         0.001 100; % (11) kappa_el
         0.001 100; % (12) D_el
         ]; 

        lb = 0.001*bounds(:,1); %바운더리 설정
        ub = 30*bounds(:,2);

        factors_ini = ones(size(bounds,1),1); %
        SOC = soc_vec(i);
        soc = SOC*0.01; % [1]

        if soc_integ == 1 && i == 1
           [f_data, z_data_integ_sep, z_model_integ_sep, z_model_dist_integ_sep, paras_integ ,paras_integ_dist] = BSL_func_3E_calc_el_multi_soc_2nd_anode(path_folder, path_file, multi_soc_range, cell_type, anode_type);
           z_data_integ = z_data_integ_sep(:,1:end/2) + z_data_integ_sep(:,end/2+1:end);
           z_model_integ = z_model_integ_sep(:,1:end/2) + z_model_integ_sep(:,end/2+1:end);
           z_model_dist_integ = z_model_dist_integ_sep(:,1:end/2) + z_model_dist_integ_sep(:,end/2+1:end);
            for j = 1:length(multi_soc_range)
                SOC = multi_soc_range(j);
                figure(j)
                t = tiledlayout(1,2,"TileSpacing","compact",'Padding','compact');
                nexttile
                plot(z_data_integ(:,2*j-1),-z_data_integ(:,2*j),'ok','linewidth',1); hold on
                plot(z_model_integ(:,2*j-1),-z_model_integ(:,2*j),'or','linewidth',1)
                plot(z_model_dist_integ(:,2*j-1),-z_model_dist_integ(:,2*j),'o','linewidth',1,'Color',[0.3010 0.7450 0.9330])
                
                legend('Exp Data','P2D',['P2D' ' + ' dist])
                axis_limit = 1.1*max(max(abs(z_data_integ(:,2*j-1:2*j))));
                    set(gca,'Box','on',... %Axis Properties: BOX   
                    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
                    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
                    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
                
                    grid on;
                    title([cell_type ' soc' num2str(SOC) ' P2D'         ' + ' dist])
                    xlabel('Z_{re} [Ohm]')
                    ylabel('-Z_{im} [Ohm]')
                hold off
                
                
                % Zoom-in semicircle
                nexttile
                plot(z_data_integ(:,2*j-1),-z_data_integ(:,2*j),'ok','linewidth',1); hold on
                plot(z_model_integ(:,2*j-1),-z_model_integ(:,2*j),'or','linewidth',1)
                plot(z_model_dist_integ(:,2*j-1),-z_model_dist_integ(:,2*j),'o','linewidth',1,'Color',[0.3010 0.7450 0.9330])
                legend('Exp Data','P2D',['P2D' ' + ' dist])
                 f_zoom_lb = 10; %[Hz] 
                    idx_zoom = f_data>f_zoom_lb;
                    axis_limit = 1.1*max(max(abs(z_data_integ(idx_zoom,:))));
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
                
                if save_check == 1
                    fig_name = sprintf([anode_type ' ' cell_type ' soc_%d' ' P2D' ' + ' dist ' .fig'],SOC); %피규어 저장
                    savefig(fullfile(save_path,fig_name));
                end

                  figure()
                t = tiledlayout(2,2,"TileSpacing","compact",'Padding','compact');
                % anode
                nexttile
                plot(z_data_integ(:,1),-(z_data_integ(:,2)),'ok','linewidth',1); hold on
                plot(z_model_integ(:,1),-z_model_integ(:,2),'or','linewidth',1)
                plot(z_model_dist_integ(:,1),-z_model_dist_integ(:,2),'o','linewidth',1,'Color',[0.3010 0.7450 0.9330])
                
                legend('Exp Data','P2D',['P2D' ' + ' dist ' anode'])
                axis_limit = 1.1*max(max(abs(z_data_integ)));
                set(gca,'Box','on',... %Axis Properties: BOX   
                'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
                'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
                'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
            
                grid on;
                title([cell_type ' anode' ' soc' num2str(SOC) ' P2D' ' + ' dist])
                xlabel('Z_{re} [Ohm]')
                ylabel('-Z_{im} [Ohm]')
                hold off
    
                % Zoom-in semicircle
                nexttile
                plot(z_data_integ(:,1),-(z_data_integ(:,2)),'ok','linewidth',1); hold on
                plot(z_model_integ(:,1),-z_model_integ(:,2),'or','linewidth',1)
                plot(z_model_dist_integ(:,1),-z_model_dist_integ(:,2),'o','linewidth',1,'Color',[0.3010 0.7450 0.9330])
                legend('Exp Data','P2D',['P2D' ' + ' dist ' anode'])
                
                f_zoom_lb = 10; %[Hz] 
                idx_zoom = f_data>f_zoom_lb;
                axis_limit = 1.1*max(max(abs(z_data_integ(idx_zoom,:))));
                set(gca,'Box','on',... %Axis Properties: BOX   
                'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
                'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
                'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
             
                grid on;
                title([cell_type ' anode' ' soc' num2str(SOC) ' P2D' ' + ' dist ' Zoom-in'])
                xlabel('Z_{re} [Ohm]')
                ylabel('-Z_{im} [Ohm]')
                hold off;
                
                %cathode
                nexttile
                plot(z_data_integ(:,3),-(z_data_integ(:,4)),'ok','linewidth',1); hold on
                plot(z_model_integ(:,3),-z_model_integ(:,4),'or','linewidth',1)
                plot(z_model_dist_integ(:,3),-z_model_dist_integ(:,4),'o','linewidth',1,'Color',[0.9330 0.7450 0.3010])
                
                legend('Exp Data','P2D',['P2D' ' + ' dist ' cathode'])
                axis_limit = 1.1*max(max(abs(z_data_integ)));
                set(gca,'Box','on',... %Axis Properties: BOX   
                'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
                'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
                'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
            
                grid on;
                title([cell_type ' cathode' ' soc' num2str(SOC) ' P2D' ' + ' dist])
                xlabel('Z_{re} [Ohm]')
                ylabel('-Z_{im} [Ohm]')
                hold off
    
                % Zoom-in semicircle
                nexttile
                plot(z_data_integ(:,3),-(z_data_integ(:,4)),'ok','linewidth',1); hold on
                plot(z_model_integ(:,3),-z_model_integ(:,4),'or','linewidth',1)
                plot(z_model_dist_integ(:,3),-z_model_dist_integ(:,4),'o','linewidth',1,'Color',[0.9330 0.7450 0.3010])
                legend('Exp Data','P2D',['P2D' ' + ' dist ' cathode'])
                
                f_zoom_lb = 10; %[Hz] 
                idx_zoom = f_data>f_zoom_lb;
                axis_limit = 1.1*max(max(abs(z_data_integ(idx_zoom,:))));
                set(gca,'Box','on',... %Axis Properties: BOX   
                'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
                'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
                'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
             
                grid on;
                title([cell_type ' cathode' ' soc' num2str(SOC) ' P2D' ' + ' dist ' Zoom-in'])
                xlabel('Z_{re} [Ohm]')
                ylabel('-Z_{im} [Ohm]')
                hold off;
                set(gcf,'position',[100 100 1000 1000])
    
    
    
                if save_check == 1
                fig_name = sprintf([anode_type ' ' cell_type ' half' ' soc_%d' ' P2D' ' + ' dist ' .fig'],SOC); %피규어 저장
                savefig(fullfile(save_path,fig_name));
                end
            end 
                
        
                if save_check == 1
                    paras_integ(2:end+1,:) = paras_integ;
                    paras_integ(1,:) = str2double(split(num2str(multi_soc_range)));
                    save([save_path filesep anode_type ' ' cell_type  ' P2D' ' + ' dist ' para_sum_integ '],'paras_integ')
        
                    paras_integ_dist(2:end+1,:) = paras_integ_dist;
                    paras_integ_dist(1,:) = str2double(split(num2str(multi_soc_range)));
                    save([save_path filesep anode_type ' ' cell_type  ' P2D' ' + ' dist ' para_sum_integ_dist '],'paras_integ_dist')
        
                    fprintf('Data saved successfully.\n');
                else
                    fprintf('Data not saved.\n')
                end
        elseif soc_integ == 1 && i > 1
               Useless(i) = i;
        elseif soc_integ == 0
        %% Structure Explained
    
            % Data and model func
            % for n - frequency points
            % w_data (n x 1) vector
            % [ Z_n_re Z_n_im Z_p_re Z_p_im] (n x 4) matrix
        
            % load EIS data
            name_file = sprintf(path_file, anode_type, SOC); 
            
            data = readmatrix([path_folder filesep name_file]);
            
            data_n = [data(:,1) -data(:,4) data(:,5)];
            data_p = [data(:,1) data(:,6) -data(:,7)];

            f_data = [data_n(:,1) data_p(:,1)]; %프리퀀시 데이터 통합
                if any(f_data(:,1)~=f_data(:,2))
                    error('freq range of EIS_p and EIS_n are different')
                end
        
             z_data = [data_n(:,2) data_n(:,3) data_p(:,2) data_p(:,3)]; %Z_real, imag 데이터 통합.
                   % [Z_n_re     Z_n_im      Z_p_re      Z_p_im]    (n x 4) matrix
             z_data_sum = [z_data(:,1) + z_data(:,3) z_data(:,2) + z_data(:,4)];
              
             figure(1)
             plot(data_n(:,2),-data_n(:,3),'ro'); hold on; grid on %음극 초도 플랏
             plot(data_p(:,2),-data_p(:,3),'bo'); %양극 초도 플랏
            
             axis_limit = 1.1*max(max(abs(z_data_sum)));%축 제한 값 설정 (최대값 1.1배로 설정)
             set(gca,'Box','on',... %Axis Properties: BOX   
             'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
             'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
             'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
             xlabel('Z_{re} [Ohm]')
             ylabel('-Z_{im} [Ohm]')
        
             % trim the high-frequency inductance part (positive imag data)
             ind_keep = data_n(:,3) <= 0 & data_p(:,3) <= 0; %0 이상인 imag data 데이터 자름. (0미만인 부분만 index 지정) 
             f_data = f_data(ind_keep); %프리퀀시도 트림
             z_data = z_data(ind_keep,:); %z_data 트림

             figure(1)
             plot(z_data(:,1),-z_data(:,2),'ro','MarkerFaceColor','red'); hold on;  %트림한 데이터 플랏 음극, 양극
             plot(z_data(:,3),-z_data(:,4),'bo','MarkerFaceColor','blue');
             hold off;
             legend('Anode','Cathode');
             grid on;
 
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
           weighted_model = @(factors,f_data)BSL_func_EISmodel_V1_3E_2nd_anode(f_data, factors,soc,T,type_acf,type_anode,type_cell).*weight_matrix; %모델 예측값에 가중치 곱하여 모델 설정
           %model = @(factors,f_data)BSL_func_EISmodel_V1_3E(f_data, factors,soc,T,type_acf);
           weighted_data = z_data.*weight_matrix; %가중치 곱해진 데이터 설정
        %   fitting
           tic;
           [factors_hat, resnorm,residual,~,~,~,jacobian_hat] ...
                = lsqcurvefit(weighted_model,factors_ini,...
                              f_data,weighted_data, lb, ub, options); %factor_ini 최적화. weighted model 및 data 사용. --> 초기점 설정 후 EIS_model 에서 계산해 파라미터 얻음
           toc;


        %% Plot Results
            [z_model1, paras1] = BSL_func_EISmodel_V1_3E_2nd_anode(f_data,factors_hat,soc,T,type_acf,type_anode,type_cell);
            
            z_data_simul1 = [z_model1(:,1)+z_model1(:,3) z_model1(:,2)+z_model1(:,4)];

            % Nyquist Plot
            % line colors
            cmat_jet = jet(16); %컬러맵 배열 생성
            cmat_n = [cmat_jet(16,:);cmat_jet(14,:);cmat_jet(11,:)]; %애노드, 캐소드 각각 색 지정
            cmat_p = [cmat_jet(1,:);cmat_jet(3,:);cmat_jet(6,:)];

            figure(2)
            t = tiledlayout(1,2); %2x2 피규어 생성, set handle at t
            nexttile
            hold on
            plot(z_data(:,1)+z_data(:,3),-(z_data(:,2)+z_data(:,4)),'o','linewidth',1,'color',cmat_n(1,:));  
            % plot(z_data_sum0(:,1),-z_data_sum0(:,2),'o','linewidth',1,'color',cmat_n(3,:))
            plot(z_data_simul1(:,1),-z_data_simul1(:,2),'o','linewidth',1,'color',cmat_n(2,:)) 
            legend('Exp Data','3E Simul')
            %legend('Exp Data','Model Fit')
            daspect ([1 1 2])
            
            
            axis_limit = 1.1*max(max(abs(z_data_sum)));
            set(gca,'Box','on',... %Axis Properties: BOX   
            'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
            'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
            'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
            grid on;
            hold off
            title(sprintf('%s 3E Simul plot at soc%d',name_type,SOC))
            xlabel('Z_{re} [Ohm]')
            ylabel('-Z_{im} [Ohm]')

            %Zoom-in semicircle
            nexttile
            hold on
            plot(z_data(:,1)+z_data(:,3),-(z_data(:,2)+z_data(:,4)),'o','linewidth',1,'color',cmat_n(1,:));  
            % plot(z_data_sum0(:,1),-z_data_sum0(:,2),'o','linewidth',1,'color',cmat_n(3,:))
            plot(z_data_simul1(:,1),-z_data_simul1(:,2),'o','linewidth',1,'color',cmat_n(2,:)) 
            legend('Exp Data','3E Simul')
            %legend('Exp Data','Model Fit')
            daspect ([1 1 2])

            f_zoom_lb = 10; %[Hz] 
            idx_zoom = f_data>f_zoom_lb;
            axis_limit = 1.1*max(max(abs(z_data_sum(idx_zoom,:))));
            set(gca,'Box','on',... %Axis Properties: BOX   
            'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
            'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
            'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
            grid on;
            hold off
            title(sprintf('%s 3E Simul plot at soc%d',name_type,SOC))
            xlabel('Z_{re} [Ohm]')
            ylabel('-Z_{im} [Ohm]')

            t.TileSpacing = 'compact';
            t.Padding = 'compact';
            set(gcf,'position',[100 100 1200 600])

            pause(0.1);
        %% Fitting Improvement by Distributed models.
            std_ini = 0.5;    
            factors_ini = [factors_hat' std_ini std_ini std_ini std_ini]; % drta,c / ddt,ac 분포 각각 고려
        
            ub = factors_ini*100;
            lb = factors_ini*0.01;

         %% [dist] Call EIS model
            weighted_model_dist = @(factors,f_data)BSL_func_EISmodel_V_3E_Dist_integrated_2nd_anode(f_data, factors,soc,T,type_acf,type_dist,type_anode,type_cell)...
               .*weight_matrix;
            weighted_data = z_data.*weight_matrix;
        
            tic;
            [factors_hat_dist, resnorm,residual,~,~,~,jacobian_hat] ...
                = lsqcurvefit(weighted_model_dist,factors_ini,...
                              f_data,weighted_data, lb, ub, options_dist);
            toc;

        %% [dist] fitting result
            % Nyquist Plot
            [z_model2, paras2] = BSL_func_EISmodel_V_3E_Dist_integrated_2nd_anode(f_data,factors_hat_dist,soc,T,type_acf,type_dist,type_anode,type_cell);
            
             z_data_simul2 = [z_model2(:,1)+z_model2(:,3) z_model2(:,2)+z_model2(:,4)];

            figure(3)
            t = tiledlayout(1,2,"TileSpacing","compact",'Padding','compact');
            nexttile
            plot(z_data(:,1)+z_data(:,3),-(z_data(:,2)+z_data(:,4)),'ok','linewidth',1); hold on
            plot(z_data_simul1(:,1),-z_data_simul1(:,2),'or','linewidth',1)
            plot(z_data_simul2(:,1),-z_data_simul2(:,2),'o','linewidth',1,'Color',[0.3010 0.7450 0.9330])
            
            legend('Exp Data','P2D',['P2D' ' + ' dist])
            axis_limit = 1.1*max(max(abs(z_data_sum)));
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
            plot(z_data(:,1)+z_data(:,3),-(z_data(:,2)+z_data(:,4)),'ok','linewidth',1); hold on
            plot(z_data_simul1(:,1),-z_data_simul1(:,2),'or','linewidth',1)
            plot(z_data_simul2(:,1),-z_data_simul2(:,2),'o','linewidth',1,'Color',[0.3010 0.7450 0.9330])
            legend('Exp Data','P2D',['P2D' ' + ' dist])
            
            f_zoom_lb = 10; %[Hz] 
            idx_zoom = f_data>f_zoom_lb;
            axis_limit = 1.1*max(max(abs(z_data_sum(idx_zoom,:))));
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


            if save_check == 1
            fig_name = sprintf([anode_type ' ' cell_type ' soc_%d' ' P2D' ' + ' dist ' .fig'],SOC); %피규어 저장
            savefig(fullfile(save_path,fig_name));
            end
            
            figure(4)
            t = tiledlayout(2,2,"TileSpacing","compact",'Padding','compact');
            % anode
            nexttile
            plot(z_data(:,1),-(z_data(:,2)),'ok','linewidth',1); hold on
            plot(z_model1(:,1),-z_model1(:,2),'or','linewidth',1)
            plot(z_model2(:,1),-z_model2(:,2),'o','linewidth',1,'Color',[0.3010 0.7450 0.9330])
            
            legend('Exp Data','P2D',['P2D' ' + ' dist ' anode'])
            axis_limit = 1.1*max(max(abs(z_data)));
            set(gca,'Box','on',... %Axis Properties: BOX   
            'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
            'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
            'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
        
            grid on;
            title([cell_type ' anode' ' soc' num2str(SOC) ' P2D' ' + ' dist])
            xlabel('Z_{re} [Ohm]')
            ylabel('-Z_{im} [Ohm]')
            hold off

            % Zoom-in semicircle
            nexttile
            plot(z_data(:,1),-(z_data(:,2)),'ok','linewidth',1); hold on
            plot(z_model1(:,1),-z_model1(:,2),'or','linewidth',1)
            plot(z_model2(:,1),-z_model2(:,2),'o','linewidth',1,'Color',[0.3010 0.7450 0.9330])
            legend('Exp Data','P2D',['P2D' ' + ' dist ' anode'])
            
            f_zoom_lb = 10; %[Hz] 
            idx_zoom = f_data>f_zoom_lb;
            axis_limit = 1.1*max(max(abs(z_data(idx_zoom,:))));
            set(gca,'Box','on',... %Axis Properties: BOX   
            'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
            'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
            'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
         
            grid on;
            title([cell_type ' anode' ' soc' num2str(SOC) ' P2D' ' + ' dist ' Zoom-in'])
            xlabel('Z_{re} [Ohm]')
            ylabel('-Z_{im} [Ohm]')
            hold off;
            
            %cathode
            nexttile
            plot(z_data(:,3),-(z_data(:,4)),'ok','linewidth',1); hold on
            plot(z_model1(:,3),-z_model1(:,4),'or','linewidth',1)
            plot(z_model2(:,3),-z_model2(:,4),'o','linewidth',1,'Color',[0.9330 0.7450 0.3010])
            
            legend('Exp Data','P2D',['P2D' ' + ' dist ' cathode'])
            axis_limit = 1.1*max(max(abs(z_data)));
            set(gca,'Box','on',... %Axis Properties: BOX   
            'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
            'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
            'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
        
            grid on;
            title([cell_type ' cathode' ' soc' num2str(SOC) ' P2D' ' + ' dist])
            xlabel('Z_{re} [Ohm]')
            ylabel('-Z_{im} [Ohm]')
            hold off

            % Zoom-in semicircle
            nexttile
            plot(z_data(:,3),-(z_data(:,4)),'ok','linewidth',1); hold on
            plot(z_model1(:,3),-z_model1(:,4),'or','linewidth',1)
            plot(z_model2(:,3),-z_model2(:,4),'o','linewidth',1,'Color',[0.9330 0.7450 0.3010])
            legend('Exp Data','P2D',['P2D' ' + ' dist ' cathode'])
            
            f_zoom_lb = 10; %[Hz] 
            idx_zoom = f_data>f_zoom_lb;
            axis_limit = 1.1*max(max(abs(z_data(idx_zoom,:))));
            set(gca,'Box','on',... %Axis Properties: BOX   
            'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
            'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
            'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
         
            grid on;
            title([cell_type ' cathode' ' soc' num2str(SOC) ' P2D' ' + ' dist ' Zoom-in'])
            xlabel('Z_{re} [Ohm]')
            ylabel('-Z_{im} [Ohm]')
            hold off;
            set(gcf,'position',[100 100 1000 1000])


            para_sum(1,1+SOC/10) = SOC; %파라미터 하나로 뭉치기
            para_sum(2:13,1+SOC/10) = paras1;
            para_sum(15:30,1+SOC/10) = paras2;
            para_sum(32,1+SOC/10) = resnorm; 

            if save_check == 1
            fig_name = sprintf([anode_type ' ' cell_type ' half' ' soc_%d' ' P2D' ' + ' dist ' .fig'],SOC); %피규어 저장
            savefig(fullfile(save_path,fig_name));
            end
            
    end
    end
            if save_check == 1
                save([save_path filesep anode_type ' ' cell_type  ' P2D' ' + ' dist ' para_sum '],'para_sum')
                fprintf('Data saved successfully.\n');
            else
                fprintf('Data not saved.\n')
            end
end

% end