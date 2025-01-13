% This code is for plotting partial frequency data of EIS
% plot the datas from soc0 to 100

clc, clear, close all;

% EIS data path
% path_folder = 'G:\공유 드라이브\BSL-Data\LGES\LG raw data\12_6cm2_soc10_EIS # Sample 1';%불러올 데이터 폴더 경로 지정 
% path_file = 'PEIS_C09_%s_cycle_soc%d.csv'; % %s: anode or cathode, %d: soc percent

% path_folder = 'G:\공유 드라이브\BSL-Data\LGES\Modeling data_2차 음극 2종\processed_data';%Blend
% path_file = '인조_천연 Blending음극_#1_soc%d_%s.csv'; % %s: anode or cathode, %d: soc percent

path_folder = 'G:\공유 드라이브\BSL-Data\LGES\Modeling data_2차 음극 2종\processed_data';%Natural
path_file = '천역흑연100_음극_#2_soc%d_%s.csv'; % %s: anode or cathode, %d: soc percent

% Partial freq setting (total 150 freq, now use 107 - 110) (blend: total
% 136, now use 109 - 112) (natural: total 136 109 - 112) 
freq_start = 109;
freq_end = 112;


% load datas
soc_range = 0:10:100;

% load whole range of datas
for i = 1:length(soc_range)
    % name_filea = sprintf(path_file,'Anode',soc_range(i)); 
    % name_filec = sprintf(path_file,'Cathode',soc_range(i));

    name_filea = sprintf(path_file,soc_range(i),'Anode'); %for anode 2 sets
    name_filec = sprintf(path_file,soc_range(i),'Cathode');

    dataa_raw(:,3*i-2:3*i) = readmatrix([path_folder filesep name_filea]);
    dataa(:,3*i-2:3*i) = dataa_raw(dataa_raw(:,3) <= 0 ,3*i-2:3*i);

    datac_raw(:,3*i-2:3*i) = readmatrix([path_folder filesep name_filec]);
    datac(:,3*i-2:3*i) = datac_raw(datac_raw(:,3) <= 0 ,3*i-2:3*i);
end

%% Bode plot
%calculation
figure()
t = tiledlayout(4,3,"TileSpacing","tight","Padding","tight");
for i = 1:length(soc_range)
    nexttile
    magnitudea = sqrt(dataa(:,3*i-1).^2 + dataa(:,3*i).^2); %|z|
    phase = atan(dataa(:,3*i)./dataa(:,3*i-1)).*180/pi; % degree phase
    
    magnitudec = sqrt(datac(:,11).^2 + datac(:,12).^2); %|z|
    hold on;
    % hold on;
    % yyaxis left
    plot(log10(dataa(:,3*i-2)), log10(magnitudea), 'bo-')
    plot(log10(datac(:,3*i-2)), log10(magnitudec), 'ro-')
    % yyaxis right
    % plot(log10(dataa(:,1)), phase, 'ro-')
    hold off;
    
    title(['SOC' num2str(soc_range(i))])
    legend('anode','cathode')
    set(gca, 'PlotBoxAspectRatio',[2.5 1 1],... 
        'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    box on;

end 
set(gcf,'Position',[100 100 1200 800])
%% Calc EIS Del bump magnitude
Dela_bump = zeros(1,length(soc_range));
Delc_bump = zeros(1,length(soc_range));

for k = 1:length(soc_range)
    Dela_bump(k) = sqrt((dataa(freq_end,3*k-1)-dataa(freq_start,3*k-1))^2 + (dataa(freq_end,3*k)-dataa(freq_start,3*k))^2); %Real & Imag in sequence
    Delc_bump(k) = sqrt((datac(freq_end,3*k-1)-datac(freq_start,3*k-1))^2 + (datac(freq_end,3*k)-datac(freq_start,3*k))^2); %Real & Imag in sequence
end