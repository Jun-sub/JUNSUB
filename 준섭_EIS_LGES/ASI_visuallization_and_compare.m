clear, clc, close all; 

%% Configuration

phi = [14 16] %Radius of Coin cell cathode & anode, [mm]

Ac = [pi*(phi(1)*10^-1/2)^2 pi*(phi(2)*10^-1/2)^2]; %Area, [cm^2]
area = [Ac(1), 12.6]; %[cm^2] 4.496 Cathode, 4.693 Anode, 1: coincell, 2: Pouch, 12.6 Cathode, 13.33 Anode 
%area = [Ac(2), 13.33]; % Anode

SOC = 4; %SOC 0 = 1, 20 = 2, 40 = 3, 60 = 4, 80 = 5, 100 = 6
SOC_num = 60; %for naming_real SOC value

% data load
folder_1 = 'G:\공유 드라이브\BSL-Data\LGES\2차 실험\RPT'; %BSL data
file_1 = 'CHC_C_20_RPT 1_01_01_EIS.txt'; % BSL data

folder_2 = 'G:\공유 드라이브\BSL-Data\LGES\LG raw data\12_6cm2_soc10_EIS # Sample 1'; %LGES data
file_2 = 'PEIS_C09_cathode_cycle_soc60.csv'; % LGES data

data_1 = readtable([folder_1 filesep file_1], 'NumHeaderLines',13,'ReadVariableNames',0); %BSL data
data_2 = readmatrix([folder_2 filesep file_2]); %LGES data

% BSL_data set
time_cell_1 = data_1.Var2;  
z_re_1 = data_1.Var12; 
z_im_1 = data_1.Var13;
cycle = data_1.Var5;

z_data = [z_re_1 z_im_1];
% LGES_data set
freq_2 = data_2(:,1);   
z_re_2 = data_2(:,2); 
z_im_2 = data_2(:,3);

idx = (cycle == SOC);

%pre-plot
figure()
hold on;

plot(z_re_1(idx), -z_im_1(idx),'ko');
plot(z_re_2, -z_im_2,'bo');

title(['EIS', ' SOC ', num2str(SOC_num)]);
xlabel ("Z'");
ylabel ("z''");
grid on;
legend ('Experimental coin cell' , 'LGES pouch cell');
hold off;

axis_limit = 1.2*max(max(abs(z_data(idx))));
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
% ASI plot
figure()
hold on;

plot(area(1)*z_re_1(idx), area(1)*-z_im_1(idx),'ro');
plot(area(2)*z_re_2, area(2)*-z_im_2,'bo');
hold off;

title(['ASI', ' SOC ', num2str(SOC_num)]);
xlabel ("Z'");
ylabel ("z''");
grid on;
legend('BSL data', 'LGES');

axis_limit = area(1)*1.2*max(max(abs(z_data(idx))));
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])