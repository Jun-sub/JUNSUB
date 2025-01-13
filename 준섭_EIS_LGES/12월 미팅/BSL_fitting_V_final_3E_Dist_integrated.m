%% 
% This function 
% (1) finds the best-fit parameters for a halfcells calls a EIS model function
% (2) improve fitting result by addapting additional distribution model
% V2 is only for 3E_Simul & multi-soc fitting
%Configuration 설정 후 최적화 진행

%결과는 para_sum으로 저장, (1) SOC P2D model 사용 (2) R_itsc (3) i0 (4) C_dl (5) Ds (6) kappa_el (7) D_el (8) Av
%P2D + Dist model 사용 (11) R_itsc (12) i0 (13) C_dl (14) Ds (15) kappa_el
%(16) D_el (17) Av (18) sigma
%---------------------------------------------------------------------------------------%



clear; clc; close all

%% Configurations
% multi-soc range 
    multi_start = 20;
    multi_end = 80;
    multi_soc_interval = 10;

% Fitting configuration
    type_weight = 1; % 0 for absoulte error, 1 for relative error

    type_anode = 2; % 0 for base, 1 for natural 2 for blend
    type_dist = 2; % 0 for DRT, 1 for DDT, 2 for integrated (DRT + DDT)

    num_iter = 100; %P2D optimization max iter
    num_iter_dist = 0; % Dist optimization max iter

% Temperature
    T = 298.15; %[K]
 
save_path = 'G:\공유 드라이브\Battery Software Lab\Projects\LGES 2023\발표 및 공유자료\최종 보고서\Kel, Del, AV_linear, y_shift separation\blending_2';
save_check = 1; % 0: don't save, 1: save
% 주의: 동일한 폴더에 동일한 type_acf, type_dist 사용시 기존 파일 삭제 후 저장됨

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
    type_acf = 4; % for 3E_Simul
    multi_soc_range = multi_start:multi_soc_interval:multi_end; % for soc integrate


       if type_dist == 0
            dist = 'soc integrated + DRT';
       elseif type_dist == 1 
            dist = 'soc integrated + DDT';
       elseif type_dist == 2
            dist = 'soc integrated + DRT + DDT';
       end


       cell_type = '3E simul';

    disp([cell_type ' ' dist]);
    fprintf('\n')
    disp(['Fitting SOC: ' num2str(multi_soc_range)])

%% 3E_Simul script
    % Optimization options
    options= optimset('Display','iter','MaxIter',num_iter,'MaxFunEvals',1e6,...
        'TolFun',1e-8,'TolX',1e-8,'FinDiffType','central');

    options_dist= optimset('Display','iter','MaxIter',num_iter_dist,'MaxFunEvals',1e6,...
        'TolFun',1e-8,'TolX',1e-8,'FinDiffType','central');

 
        % Parameters 
        %   1         2     3     4      5     6         7      8     9   10    11     12 
        % R_itsc_n, i0_n, Cdl_n, Ds_n, Av_n, R_itsc_p, i0_p, Cdl_p, Ds_p, Av_p, k_el, D_el
        
        bounds = [...
             0.5 2 % (1) R_n
             0.1 50; % (2) i0_n
             0.1 50; % (3) C_dl_n
             0.01 100; % (4) Ds_n
             0.02 10; % (5) Av_n
             0.5 2;  % (6) R_p
             0.1 50; % (7) i0_p
             0.1 50; % (8) C_dl_p
             0.01 100; % (9) Ds_p
             0.1 10; % (10) Av_p
             0.1 10; % (11) kappa_el
             0.001 100; % (12) D_el
             ]; 

       
        lb = 0.1*bounds(:,1); %바운더리 설정
        ub = 30*bounds(:,2);

        factors_ini = ones(1,size(bounds,1)); %

        % Start multi-soc fitting

           [f_data, z_data_integ_sep, z_model_integ_sep, z_model_dist_integ_sep, paras_integ ,paras_integ_dist] = BSL_func_3E_calc_el_multi_soc_V_final(path_folder, path_file, multi_soc_range, cell_type, []);
           
           z_data_integ = z_data_integ_sep(:,1:end/2) + z_data_integ_sep(:,end/2+1:end);
           z_model_integ = z_model_integ_sep(:,1:end/2) + z_model_integ_sep(:,end/2+1:end);
           z_model_integ_dist = z_model_dist_integ_sep(:,1:end/2) + z_model_dist_integ_sep(:,end/2+1:end);

           figure(length(multi_soc_range)*3) % for preventing of figure number overlap with main figure / Anode
           t = tiledlayout(ceil(sqrt(length(multi_soc_range))),ceil(sqrt(length(multi_soc_range))),"TileSpacing","loose","Padding","loose");
            
           for i = 1:length(multi_soc_range)
                nexttile 
                hold on;
                plot(z_data_integ_sep(:,2*i-1),-z_data_integ_sep(:,2*i),'ok','linewidth',1,'MarkerSize',4)
                plot(z_model_integ_sep(:,2*i-1),-z_model_integ_sep(:,2*i),'ob-','linewidth',1,'MarkerSize',6)
                plot(z_model_dist_integ_sep(:,2*i-1),-z_model_dist_integ_sep(:,2*i),'o-','linewidth',1,'Color',[0.3010 0.7450 0.9330],'MarkerSize',6)
                
                axis_limit = 1.1*max(max(abs(z_data_integ_sep(:,2*i-1:2*i))));
                legend(['z data' ' soc ' num2str(multi_soc_range(i))],['z model' ' soc ' num2str(multi_soc_range(i))],['z model + DRT + DDT' ' soc ' num2str(multi_soc_range(i))])
                title('Anode')
                set(gca,'Box','on',... %Axis Properties: BOX   
                'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
                'FontUnits','points','FontSize',10,'FontName','Times New Roman','XLim',[0 axis_limit],'Ylim',[0 axis_limit])
                grid on;
                xlabel('Z_{re} [Ohm]')
                ylabel('-Z_{im} [Ohm]')
                hold off;
            end

            set(gcf,'Position',[100 100 1200 1200]);
            figure(length(multi_soc_range)*4) % for preventing of figure number overlap with main figure / Cathode
             t = tiledlayout(ceil(sqrt(length(multi_soc_range))),ceil(sqrt(length(multi_soc_range))),"TileSpacing","loose","Padding","loose");
                for i = 1:length(multi_soc_range)
                    nexttile 
                    hold on;
                    plot(z_data_integ_sep(:,2*(i+length(multi_soc_range))-1),-z_data_integ_sep(:,2*(i+length(multi_soc_range))),'ok','linewidth',1,'MarkerSize',4)
                    plot(z_model_integ_sep(:,2*(i+length(multi_soc_range))-1),-z_model_integ_sep(:,2*(i+length(multi_soc_range))),'or-','linewidth',1,'MarkerSize',6)
                    plot(z_model_dist_integ_sep(:,2*(i+length(multi_soc_range))-1),-z_model_dist_integ_sep(:,2*(i+length(multi_soc_range))),'o-','linewidth',1,'Color',[0.9330 0.7450 0.3010],'MarkerSize',6)

                    axis_limit = 1.1*max(max(abs(z_data_integ_sep(:,2*(i+length(multi_soc_range))-1:2*(i+length(multi_soc_range))))));
                    legend(['z data' ' soc ' num2str(multi_soc_range(i))],['z model' ' soc ' num2str(multi_soc_range(i))],['z model + DRT + DDT' ' soc ' num2str(multi_soc_range(i))])
                    title('Cathode')
                    set(gca,'Box','on',... %Axis Properties: BOX   
                    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
                    'FontUnits','points','FontSize',10,'FontName','Times New Roman','XLim',[0 axis_limit],'Ylim',[0 axis_limit])
                    grid on;
                    xlabel('Z_{re} [Ohm]')
                    ylabel('-Z_{im} [Ohm]')
                    hold off;
                end

           set(gcf,'Position',[1300 100 1200 1200]);
           pause(0.1); %for exhibit figure

           if save_check == 1
              figure(length(multi_soc_range)*3)
              savefig(fullfile(save_path, sprintf(['Anode' ' P2D' ' + ' dist ' integ' ' .fig'])));

              figure(length(multi_soc_range)*4)
              savefig(fullfile(save_path, sprintf(['Cathode' ' P2D' ' + ' dist ' integ' '.fig'])));

              save(fullfile(save_path, sprintf([cell_type ' P2D' ' + ' dist ' integ workspace' '.mat'])));
           end 

           
           % for plotting each soc figure

            for j = 1:length(multi_soc_range)
                SOC = multi_soc_range(j);
                figure(j) % for 3E_Simul
                t = tiledlayout(1,2,"TileSpacing","compact",'Padding','compact');
                nexttile
                plot(z_data_integ(:,2*j-1),-z_data_integ(:,2*j),'ok','linewidth',1); hold on
                plot(z_model_integ(:,2*j-1),-z_model_integ(:,2*j),'or','linewidth',1)
                plot(z_model_integ_dist(:,2*j-1),-z_model_integ_dist(:,2*j),'o','linewidth',1,'Color',[0.3010 0.7450 0.9330])

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
                plot(z_model_integ_dist(:,2*j-1),-z_model_integ_dist(:,2*j),'o','linewidth',1,'Color',[0.3010 0.7450 0.9330])
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
                    fig_name = sprintf([cell_type ' soc_%d' ' P2D' ' + ' dist ' .fig'],SOC); %피규어 저장
                    savefig(fullfile(save_path,fig_name));
                end


                figure(length(multi_soc_range)+j)
                t = tiledlayout(2,2,"TileSpacing","compact",'Padding','compact');
                % anode
                nexttile
                plot(z_data_integ_sep(:,2*j-1),-z_data_integ_sep(:,2*j),'ok','linewidth',1); hold on
                plot(z_model_integ_sep(:,2*j-1),-z_model_integ_sep(:,2*j),'ob','linewidth',1)
                plot(z_model_dist_integ_sep(:,2*j-1),-z_model_dist_integ_sep(:,2*j),'o','linewidth',1,'Color',[0.3010 0.7450 0.9330])

                legend('Exp Data','P2D',['P2D' ' + ' dist ' anode'])
                axis_limit = 1.1*max(max(abs(z_data_integ_sep(:,2*j-1:2*j))));
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
                plot(z_data_integ_sep(:,2*j-1),-z_data_integ_sep(:,2*j),'ok','linewidth',1); hold on
                plot(z_model_integ_sep(:,2*j-1),-z_model_integ_sep(:,2*j),'ob','linewidth',1)
                plot(z_model_dist_integ_sep(:,2*j-1),-z_model_dist_integ_sep(:,2*j),'o','linewidth',1,'Color',[0.3010 0.7450 0.9330])
                legend('Exp Data','P2D',['P2D' ' + ' dist ' anode'])

                f_zoom_lb = 10; %[Hz] 
                idx_zoom = f_data>f_zoom_lb;
                axis_limit = 1.1*max(max(abs(z_data_integ_sep(idx_zoom,2*j-1:2*j))));
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
                plot(z_data_integ_sep(:,2*(length(multi_soc_range)+j)-1),-z_data_integ_sep(:,2*(length(multi_soc_range)+j)),'ok','linewidth',1); hold on
                plot(z_model_integ_sep(:,2*(length(multi_soc_range)+j)-1),-z_model_integ_sep(:,2*(length(multi_soc_range)+j)),'or','linewidth',1)
                plot(z_model_dist_integ_sep(:,2*(length(multi_soc_range)+j)-1),-z_model_dist_integ_sep(:,2*(length(multi_soc_range)+j)),'o','linewidth',1,'Color',[0.9330 0.7450 0.3010])

                legend('Exp Data','P2D',['P2D' ' + ' dist ' cathode'])
                axis_limit = 1.1*max(max(abs(z_data_integ_sep(:,2*(length(multi_soc_range)+j)-1:2*(length(multi_soc_range)+j)))));
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
               plot(z_data_integ_sep(:,2*(length(multi_soc_range)+j)-1),-z_data_integ_sep(:,2*(length(multi_soc_range)+j)),'ok','linewidth',1); hold on
                plot(z_model_integ_sep(:,2*(length(multi_soc_range)+j)-1),-z_model_integ_sep(:,2*(length(multi_soc_range)+j)),'or','linewidth',1)
                plot(z_model_dist_integ_sep(:,2*(length(multi_soc_range)+j)-1),-z_model_dist_integ_sep(:,2*(length(multi_soc_range)+j)),'o','linewidth',1,'Color',[0.9330 0.7450 0.3010])
                legend('Exp Data','P2D',['P2D' ' + ' dist ' cathode'])

                f_zoom_lb = 10; %[Hz] 
                idx_zoom = f_data>f_zoom_lb;
                axis_limit = 1.1*max(max(abs(z_data_integ_sep(idx_zoom,2*(length(multi_soc_range)+j)-1:2*(length(multi_soc_range)+j)))));
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
                    fig_name = sprintf([cell_type ' half' ' soc_%d' ' P2D' ' + ' dist ' .fig'],SOC); %피규어 저장
                    savefig(fullfile(save_path,fig_name));
                end

            end

            
        
                if save_check == 1
                    paras_integ(2:end+1,:) = paras_integ;
                    paras_integ(1,:) = str2double(split(num2str(multi_soc_range)));
                    save([save_path filesep cell_type  ' P2D' ' + ' dist ' para_sum_integ '],'paras_integ')
        
                    paras_integ_dist(2:end+1,:) = paras_integ_dist;
                    paras_integ_dist(1,:) = str2double(split(num2str(multi_soc_range)));
                    save([save_path filesep cell_type  ' P2D' ' + ' dist ' para_sum_integ_dist '],'paras_integ_dist')
        
                    fprintf('Data saved successfully.\n');
                else
                    fprintf('Data not saved.\n')
                end
