% This code if for plotting multi-soc upon the change of each factors

clc, clear, close all

for n = 8
% Configuration
para_range = 0.1:0.1:0.6; %3.7;
para_var = n; %1 = R_itsc, 2 = i0, 3 = Cdl, 4 = Ds, 5 = Av, 6 = DRT_std, 7 = DDT_std, 8 = Epsla, 9 = Kel, 10 = Del
para_name = {'R_itsc' 'i0' 'Cdl' 'Ds' 'Av' 'DRT_std' 'DDT_std' 'Epsla' 'Kel' 'Del'};
%data load
load("C:\Users\admin\Documents\GitHub\JunSub\준섭_EIS_LGES\10월 미팅\Av_integ_V2\work_space2.mat");
% save_path = 'C:\Users\admin\Documents\GitHub\JunSub\준섭_EIS_LGES\10월 미팅\multi-soc plot upon paras';
% plot model


% Factors correction
% if factors_integ_hat(5,1:2*length(multi_soc_range)) == 1
%    factors_integ_hat(5,1:length(multi_soc_range)) = factors_integ_hat(3,2*length(multi_soc_range)+1);
%    factors_integ_hat(5,length(multi_soc_range)+1:2*length(multi_soc_range)) = factors_integ_hat(4,2*length(multi_soc_range)+1);
%    factors_integ_hat(3:4,2*length(multi_soc_range)+1) = 0;
% end
% 
% if factors_integ_hat_dist(5,1:2*length(multi_soc_range)) == 1
%    factors_integ_hat_dist(5,1:length(multi_soc_range)) = factors_integ_hat_dist(3,2*length(multi_soc_range)+1);
%    factors_integ_hat_dist(5,length(multi_soc_range)+1:2*length(multi_soc_range)) = factors_integ_hat_dist(4,2*length(multi_soc_range)+1);
%    factors_integ_hat_dist(3:4,2*length(multi_soc_range)+1) = 0;
% end 

% P2D+DRT+DDT multi-soc
    [z_model_dist, paras_integ_dist] = BSL_func_EISmodel_V_3E_soc_and_Dist_integrated(f_data,factors_integ_hat_dist,multi_soc_range,T,type_acf,cell_type,type_dist);
    % [z_model, paras_integ] = BSL_func_soc_integrated_model_V1_3E(f_data,factors_integ_hat,multi_soc_range,T,type_acf,cell_type); 

    factors_temp = factors_integ_hat_dist;
    factors_temp2 = factors_integ_hat;
    figure() % for preventing of figure number overlap with main figure / Anode
     sgtitle('Anode')
        t = tiledlayout(ceil(sqrt(length(multi_soc_range))),ceil(sqrt(length(multi_soc_range))),"TileSpacing","tight","Padding","tight");
            tic;
            for i = 1:length(multi_soc_range)
                nexttile 
                legend_name = {['z data' ' soc ' num2str(multi_soc_range(i))],['z model dist' ' soc ' num2str(multi_soc_range(i))]};

                hold on;
                plot(z_data_integ_sep(:,2*i-1),-z_data_integ_sep(:,2*i),'ko','LineWidth',1) 
                % plot(z_model(:,2*i-1),-z_model(:,2*i),'Color',[0.7450
                % 0.9330 0.3010 ],'LineWidth',1,'Marker','o','LineStyle','none') 
                plot(z_model_dist(:,2*i-1),-z_model_dist(:,2*i),'Color',[0.3010 0.7450 0.9330],'LineWidth',1,'Marker','o','LineStyle','none')
                

                for k = 1:length(para_range)
                    if para_var <= 4 || para_var == 6
                       factors_temp(para_var,1:2*length(multi_soc_range)) = para_range(k)*factors_integ_hat_dist(para_var,1:2*length(multi_soc_range));
                       factors_temp2(para_var,1:2*length(multi_soc_range)) = para_range(k)*factors_integ_hat(para_var,1:2*length(multi_soc_range));
                    elseif para_var == 5 && factors_integ_hat(3,2*length(multi_soc_range)+1) ~= 0
                       factors_temp(3:4,2*length(multi_soc_range)+1) = para_range(k)*factors_integ_hat_dist(3:4,2*length(multi_soc_range)+1);
                    elseif para_var == 7
                       factors_temp(para_var,1:2*length(multi_soc_range)) = para_range(k);
                    elseif para_var == 8
                       factors_temp(5,2*length(multi_soc_range)+1) = para_range(k);
                    elseif para_var == 9
                       factors_temp(1,2*length(multi_soc_range)+1) = para_range(k)*factors_integ_hat_dist(1,2*length(multi_soc_range)+1);
                       factors_temp2(1,2*length(multi_soc_range)+1) = para_range(k)*factors_integ_hat(1,2*length(multi_soc_range)+1);
                    elseif para_var == 10
                       factors_temp(2,2*length(multi_soc_range)+1) = para_range(k); %*factors_integ_hat_dist(2,2*length(multi_soc_range)+1);
                       factors_temp2(2,2*length(multi_soc_range)+1) = para_range(k)*factors_integ_hat(2,2*length(multi_soc_range)+1);
                    end 
                    [z_model_multi_dist, paras_integ_dist] = BSL_func_EISmodel_V_3E_soc_and_Dist_integrated(f_data,factors_temp,multi_soc_range,T,type_acf,cell_type,type_dist);
                    % [z_model_multi, paras_integ] = BSL_func_soc_integrated_model_V1_3E(f_data,factors_temp2,multi_soc_range,T,type_acf,cell_type); 
                    plot(z_model_multi_dist(:,2*i-1),-z_model_multi_dist(:,2*i));
                    % plot(z_model_multi(:,2*i-1),-z_model_multi(:,2*i));

                    legend_name{k+2} = [para_name{para_var} ' multiply ' num2str(para_range(k))];
                    if para_var == 7 || para_var == 8 
                       legend_name{k+2} = [para_name{para_var} num2str(para_range(k))];
                    end 
                    disp(['Anode tile_' num2str(i) ' para ' num2str(k)])
                end

                legend(legend_name, 'Location','best')
                title(['Anode' ' soc ' num2str(multi_soc_range(i))])
                axis_limit = 1.1*max(max(abs(z_data_integ_sep(:,2*i-1:2*i))));
                set(gca,'Box','on',... %Axis Properties: BOX   
                'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
                'FontUnits','points','FontSize',10,'FontName','Times New Roman', 'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
                grid on;
                xlabel('Z_{re} [Ohm]')
                ylabel('-Z_{im} [Ohm]')
                hold off;

                fprintf('Process %d of %d \n\n',i,2*length(multi_soc_range))
                
            end
    
             set(gcf,'Position',[100 100 1200 1200]);

             % savefig(fullfile(save_path,['multi_soc_anode_plot_' para_name{para_var}]))


            figure() % for preventing of figure number overlap with main figure / Anode
             sgtitle('Cathode')
                t = tiledlayout(ceil(sqrt(length(multi_soc_range))),ceil(sqrt(length(multi_soc_range))),"TileSpacing","tight","Padding","tight");
            for i = 1:length(multi_soc_range)
                nexttile 
                legend_name = {['z data' ' soc ' num2str(multi_soc_range(i))],['z model dist' ' soc ' num2str(multi_soc_range(i))]};

                hold on;
                plot(z_data_integ_sep(:,2*(i+length(multi_soc_range))-1),-z_data_integ_sep(:,2*(i+length(multi_soc_range))),'ko','LineWidth',1) 
                plot(z_model_dist(:,2*(i+length(multi_soc_range))-1),-z_model_dist(:,2*(i+length(multi_soc_range))),'Color',[0.3010 0.7450 0.9330],'LineWidth',1,'Marker','o','LineStyle','none')
                
                for k = 1:length(para_range)
                    if para_var <= 4 || para_var == 6
                    factors_temp(para_var,1:2*length(multi_soc_range)) = para_range(k)*factors_integ_hat_dist(para_var,1:2*length(multi_soc_range));
                    elseif para_var == 5 && factors_integ_hat(3,2*length(multi_soc_range)+1) ~= 0
                       factors_temp(3:4,2*length(multi_soc_range)+1) = para_range(k)*factors_integ_hat_dist(3:4,2*length(multi_soc_range)+1);
                    elseif para_var == 7
                       factors_temp(para_var,1:2*length(multi_soc_range)) = para_range(k);
                    elseif para_var == 8
                       factors_temp(5,2*length(multi_soc_range)+1) = para_range(k);
                    elseif para_var == 9
                       factors_temp(1,2*length(multi_soc_range)+1) = para_range(k)*factors_integ_hat_dist(1,2*length(multi_soc_range)+1);;
                    elseif para_var == 10
                       factors_temp(2,2*length(multi_soc_range)+1) = para_range(k)*factors_integ_hat_dist(2,2*length(multi_soc_range)+1);;
                    end 

                    [z_model_multi_dist, paras_integ_dist] = BSL_func_EISmodel_V_3E_soc_and_Dist_integrated(f_data,factors_temp,multi_soc_range,T,type_acf,cell_type,type_dist);
                    plot(z_model_multi_dist(:,2*(i+length(multi_soc_range))-1),-z_model_multi_dist(:,2*(i+length(multi_soc_range))));

                    legend_name{k+2} = [para_name{para_var} ' multiply ' num2str(para_range(k))];

                    if para_var == 7 || para_var == 8 
                       legend_name{k+2} = [para_name{para_var} num2str(para_range(k))];
                    end

                    disp(['Cathode tile_' num2str(i) ' para ' num2str(k)])
                end

                legend(legend_name, 'Location','best')
                title(['Cathode' ' soc ' num2str(multi_soc_range(i))])
                axis_limit = 1.1*max(max(abs(z_data_integ_sep(:,2*(i+length(multi_soc_range))-1:2*(i+length(multi_soc_range))))));
                set(gca,'Box','on',... %Axis Properties: BOX   
                'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
                'FontUnits','points','FontSize',10,'FontName','Times New Roman', 'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
                grid on;
                xlabel('Z_{re} [Ohm]')
                ylabel('-Z_{im} [Ohm]')
                hold off;

                fprintf('Process %d of %d \n\n',length(multi_soc_range)+i,2*length(multi_soc_range))
                
            end

       
            set(gcf,'Position',[1300 100 1200 1200]);
      
            % savefig(fullfile(save_path,['multi_soc_cathode_plot_' para_name{para_var}]))
end