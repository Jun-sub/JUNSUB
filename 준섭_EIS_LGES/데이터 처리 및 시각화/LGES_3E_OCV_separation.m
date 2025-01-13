%This script is for separating raw rOCV data from 3E LGES cell.

clc, clear, close all;

folder_path = 'G:\공유 드라이브\BSL-Data\LGES\음극 2종 OCV data';
file_name = 'Natural_Gra100%_anode_ocv 추출용.txt';

data = readmatrix(fullfile(folder_path,file_name));

% data load

data_full = data(:,3);
data_anode = data(:,4);
data_cathode = data(:,5);
data_current = data(:,6);

cumQ = cumtrapz(abs(data_current));
Q = trapz(abs(data_current));

% OCV, OCP data plot
figure(1)
t = tiledlayout(4,1,"TileSpacing","compact","Padding","compact");

nexttile
plot(data_full)
legend('full')

nexttile
plot(data_anode)
legend('anode')

nexttile
plot(data_cathode)
legend('cathode')

nexttile
plot(data_current)
legend('current')

set(gcf,"Position",[200 200 800 1000])
% discharging step identification
idx_dischg = data_current ~= 0;
idx_rest = data_current == 0;

num_rest = sum(idx_rest(1:end-1) < idx_rest(2:end)); %find the number of rest step
idx_rest = find(idx_rest(1:end-1) > idx_rest(2:end)); % find the last step idx of rest (for rOCV data)

% set rOCV & rOCP data
OCV = data_full(idx_rest);
OCPa = data_anode(idx_rest);
OCPc = data_cathode(idx_rest);

soc_raw = (Q-cumQ(idx_rest))/Q;
figure(1)

nexttile(1)
hold on;
plot(idx_rest,OCV,'ko','MarkeredgeColor','black','MarkerSize',3);
hold off

nexttile(2)
hold on;
plot(idx_rest,OCPa,'ro','MarkeredgeColor','red','MarkerSize',3);
hold off

nexttile(3)
hold on;
plot(idx_rest,OCPc,'bo','MarkeredgeColor','blue','MarkerSize',3);
hold off

%plot rearrange upon soc
figure(2)
hold on;

yyaxis left
plot(soc_raw,OCPa,'ro-','MarkerSize',3)
ylabel('V')

yyaxis right
plot(soc_raw,OCPc,'bo-','MarkerSize',3)
plot(soc_raw,OCV,'ko-','MarkerSize',3)

ylabel('V')
xlabel('soc')
legend('OCPa','OCPc','OCV')
grid on;
box on;
title('SOC scale plot')
set(gca,"PlotBoxAspectRatio",[1 1 1])
set(gcf,"Position",[200 200 800 800])

%save_data
OCV(:,2) = OCV;
OCV(:,1) = soc_raw;

OCPa(:,2) = OCPa;
OCPa(:,1) = soc_raw;

OCPc(:,2) = OCPc;
OCPc(:,1) = soc_raw;

% writematrix(flip(OCV),'G:\공유 드라이브\BSL-Data\LGES\음극 2종 OCV data\추출 OCP data\Blend_OCV.csv')
% writematrix(flip(OCPa),'G:\공유 드라이브\BSL-Data\LGES\음극 2종 OCV data\추출 OCP data\Blend_OCPa.csv')
% writematrix(flip(OCPc),'G:\공유 드라이브\BSL-Data\LGES\음극 2종 OCV data\추출 OCP data\Blend_OCPc.csv')
