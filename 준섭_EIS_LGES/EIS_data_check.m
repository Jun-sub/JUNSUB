% This code is for plotting all the range of EIS datas

clc, clear, close all;

%data load
folder_path = 'G:\공유 드라이브\BSL-Data\LGES\LG raw data\12_6cm2_soc10_EIS # Sample 1';
file_name = 'PEIS_C09_%s_cycle_soc%d.csv';

% path_folder = 'G:\공유 드라이브\BSL-Data\LGES\LG raw data\12_6cm2_soc10_EIS # Sample 1';%불러올 데이터 폴더 경로 지정
% path_file = 'PEIS_C09_%s_cycle_soc%d.csv'; % %s: anode or cathode, %d: soc percent

% path_folder = 'G:\공유 드라이브\BSL-Data\LGES\Modeling data_2차 음극 2종\processed_data2';%불러올 데이터 폴더 경로 지정
% path_file = '%s_인조_천연 Blending음극_#1_soc%d.csv'; % %s: anode or cathode, %d: soc percent
% 
% 
% % for anode Natural
% path_folder = 'G:\공유 드라이브\BSL-Data\LGES\Modeling data_2차 음극 2종\processed_data2';%불러올 데이터 폴더 경로 지정
% path_file = '%s_천역흑연100_음극_#2_soc%d.csv'; % %s: anode or cathode, %d: soc percent

soc_range = 20:10:80;
type_acf = 'anode';

legend_name = cell(1,length(soc_range));
figure()
hold on;
for i = 1:length(soc_range)
    file = fullfile(folder_path,sprintf(file_name,type_acf,soc_range(i)));
    data = readmatrix(file);

    freq = data(:,1);
    z_real = data(:,2);
    z_imag = data(:,3);

    plot(z_real,-z_imag,'Marker','o')
    legend_name{i} = sprintf(['Base' '\\_' type_acf '\\_' 'soc %s'],num2str(soc_range(i)));
end
hold off
legend(legend_name)
xlabel("z'")
ylabel('z"')
box on;
grid on;
set(gca,"PlotBoxAspectRatio",[1 1 1],'FontUnits','points','FontSize',12,'FontName','TimesNewRoman','XLim',[0 6] ...
    ,'YLim',[0,6]);
set(gcf,'position',[300 300 800 800])