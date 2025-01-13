clc, clear, close all;

% Configuration
path_1 = 'G:\공유 드라이브\BSL-Data\LGES\LG raw data\12_6cm2_soc10_EIS # Sample 1'; % base
path_2 = 'G:\공유 드라이브\BSL-Data\LGES\Modeling data_2차 음극 2종\processed_data'; % Blend
path_3 = 'G:\공유 드라이브\BSL-Data\LGES\Modeling data_2차 음극 2종\processed_data'; % Natural

name_1 = 'PEIS_C09_%s_cycle_soc%d.csv'; % base
name_2 = '인조_천연 Blending음극_#1_soc%d_%s.csv'; % blend
name_3 = '천역흑연100_음극_#2_soc%d_%s.csv'; % natural

% save_path = 

type_acf = 1; % 1:anode, 2:cathode
soc = 60:10:80; %[30 50 70];

if type_acf == 1
   cell_type = 'anode';
elseif type_acf == 2
   cell_type = 'cathode';
end

% load data
for i = 1:length(soc)
    data_1(:,3*i-2:3*i) = readmatrix(fullfile(path_1,sprintf(name_1,cell_type,soc(i)))); %base data
    data_2(:,3*i-2:3*i) = readmatrix(fullfile(path_2,sprintf(name_2,soc(i),cell_type))); %blend data
    data_3(:,3*i-2:3*i) = readmatrix(fullfile(path_3,sprintf(name_3,soc(i),cell_type))); %natural data
end

% plot figure
% % base 1) each soc, 2) every DOE in 1 fig
% figure()
% t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
% 
%     for k = 1:3
%         nexttile
%         plot(data_1(:,3*k-1),-data_1(:,3*k),'Marker','o','MarkerEdgeColor','green','Color','green');
% 
%         axis_limit = 1.1*max(max(abs(data_1(:,3*k-1:3*k))));
%             set(gca,'Box','on',... %Axis Properties: BOX   
%             'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
%             'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
%             'XLim',[0 axis_limit],'Ylim',[0 axis_limit])%'XLim',[0.2 0.4],'Ylim',[0 0.2])
% 
%             title(sprintf('Base %s at SOC %s',cell_type, num2str(soc(k))))
%             grid on;
%             xlabel('Z_{re} [Ohm]')
%             ylabel('-Z_{im} [Ohm]')
%     end 
%     set(gcf,'position',[100 300 1200 400])
% 
% % blend
% figure()
% t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
% 
%     for k = 1:3
%         nexttile
%         plot(data_2(:,3*k-1),-data_2(:,3*k),'Marker','o','MarkerEdgeColor','blue','Color','blue');
% 
%         axis_limit = 1.1*max(max(abs(data_2(:,3*k-1:3*k))));
%             set(gca,'Box','on',... %Axis Properties: BOX   
%             'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
%             'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
%             'XLim',[0 axis_limit],'Ylim',[0 axis_limit])%'XLim',[0.1 0.2],'Ylim',[0 0.1])
% 
%             title(sprintf('Blend %s at SOC %s',cell_type, num2str(soc(k))))
%             grid on;
%             xlabel('Z_{re} [Ohm]')
%             ylabel('-Z_{im} [Ohm]')
%     end 
%     set(gcf,'position',[100 300 1200 400])
% 
% 
% % natural
% figure()
% t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
% 
%     for k = 1:3
%         nexttile
%         plot(data_3(:,3*k-1),-data_3(:,3*k),'Marker','o','MarkerEdgeColor','red','Color','red');
% 
%         axis_limit = 1.1*max(max(abs(data_3(:,3*k-1:3*k))));
%             set(gca,'Box','on',... %Axis Properties: BOX   
%             'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
%             'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
%             'XLim',[0 axis_limit],'Ylim',[0 axis_limit])%'XLim',[0.1 0.2],'Ylim',[0 0.1])
% 
%             title(sprintf('Natural %s at SOC %s',cell_type, num2str(soc(k))))
%             grid on;
%             xlabel('Z_{re} [Ohm]')
%             ylabel('-Z_{im} [Ohm]')
%     end 
%     set(gcf,'position',[100 300 1200 400])

% total in 1 fig
figure()
t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    
    for k = 1:3
        nexttile
        hold on;
        
        plot(data_1(:,3*k-1),-data_1(:,3*k),'Marker','o','MarkerEdgeColor','green','Color','green');
        plot(data_2(:,3*k-1),-data_2(:,3*k),'Marker','o','MarkerEdgeColor','blue','Color','blue');
        plot(data_3(:,3*k-1),-data_3(:,3*k),'Marker','o','MarkerEdgeColor','red','Color','red');
        
        axis_limit = 1.4*max(max(abs(data_1(:,3*k-1:3*k))));
            set(gca,'Box','on',... %Axis Properties: BOX   
            'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
            'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
            'XLim',[0 axis_limit],'Ylim',[0 axis_limit])%'XLim',[0.1 0.4],'Ylim',[0 0.3])
        
        title(sprintf('Total %s at SOC %s',cell_type, num2str(soc(k))))
        legend('Bese','Blend','Natural')
        grid on;
        xlabel('Z_{re} [Ohm]')
        ylabel('-Z_{im} [Ohm]')
        hold off;

        nexttile(3+k) %Zoom-in
        hold on;
        plot(data_1(:,3*k-1),-data_1(:,3*k),'Marker','o','MarkerEdgeColor','green','Color','green');
        plot(data_2(:,3*k-1),-data_2(:,3*k),'Marker','o','MarkerEdgeColor','blue','Color','blue');
        plot(data_3(:,3*k-1),-data_3(:,3*k),'Marker','o','MarkerEdgeColor','red','Color','red');
        
        axis_limit = 0.15*max(max(abs(data_1(:,3*k-1:3*k))));
            set(gca,'Box','on',... %Axis Properties: BOX   
            'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
            'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
            'XLim',[0 axis_limit],'Ylim',[0 axis_limit])%'XLim',[0.1 0.4],'Ylim',[0 0.3])
        
        title(sprintf('Total %s at SOC %s Zoon-in',cell_type, num2str(soc(k))))
        legend('Bese','Blend','Natural')
        grid on;
        xlabel('Z_{re} [Ohm]')
        ylabel('-Z_{im} [Ohm]')
        hold off;
    end 
    set(gcf,'position',[100 300 1200 800])

%Del magnitude calc
% Z1 = 0.335115; %real 1
% Z2 = 0.348797; %real 2
% 
% Z3 = 0.037189; %imag 1
% Z4 = 0.0567313; %imag 2
% 
% magnitude = sqrt((Z1 - Z2)^2 + (Z3 - Z4)^2)