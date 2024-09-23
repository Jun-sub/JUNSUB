clear; clc; close all;

% folder = 'G:\공유 드라이브\BSL-Data\LGES\2차 실험\OCP\데이터 변환';
% filename_ocpp = 'Processed_Data_CHC.mat';
% filename_ocpn = 'Processed_Data_AHC.mat';

folder = 'G:\공유 드라이브\EE6211_2024\Data\pOCV_and_pOCP';
filename_ocpn = 'AHC_(5)_OCV_C20.mat'; %parse 된 파일 사용
filename_ocpp = 'CHC_(5)_OCV_C20.mat'; 

ocpp = load([folder filesep filename_ocpp]);
ocpn = load([folder filesep filename_ocpn]);

ocpp = ocpp.OCV_golden.OCVdis; %OCPp 측정 방향 따라 OCV도 맞춤.
ocpn = ocpn.OCV_golden.OCVchg;

figure(1)
plot(ocpp(:,1),ocpp(:,2),'b-')
title('OCPp')
xlabel('Stoichiometry')
ylabel('Voltage [V]')
figure(2)
plot(ocpn(:,1),ocpn(:,2),'r-')
title('OCPn')
xlabel('Stoichiometry')
ylabel('Voltage [V]')

y0 = 0.8661%ocpp(end,1); % cathode stoi range
y1 = 0.2157%ocpp(1,1);
x0 = 0.0288%ocpn(end,1); % anode stoi range
x1 = 0.9102%ocpn(1,1);

soc = 0:0.01:1; %Interp 범위 지정

y_soc01 = y0 + (y1-y0)*soc; %Cathode, 방전시 stoi 증가하는 방향
x_soc01 = x0 + (x1-x0)*soc; %Anode, Stoi 감소하는 방향

ocpp_soc01 = interp1(ocpp(:,1),ocpp(:,2),y_soc01);
ocpn_soc01 = interp1(ocpn(:,1),ocpn(:,2),x_soc01);

ocv_soc01 = ocpp_soc01 - ocpn_soc01;

figure(3)
yyaxis left
plot(soc, ocv_soc01, 'k'); hold on; %SOC 0-1 범위 내에 OCV 표현
plot(soc, ocpp_soc01, 'b') % SOC 0-1 범위 내에 OCPp 표현
ylabel('Voltage[V]')
yyaxis right
plot(soc, ocpn_soc01, 'r');
ylabel('Voltage')
title('OCV from OCP')
xlabel('SOC')

soc = -0.2:0.01:1.2; % soc 범위 stoi 까지 보일 수 있게 확장

y_vec = y0 + (y1 - y0)*soc;
x_vec = x0 + (x1 - x0)*soc;


ocpp_vec = interp1(ocpp(:,1),ocpp(:,2),y_vec);
ocpn_vec = interp1(ocpn(:,1),ocpn(:,2),x_vec);

ocv_vec = ocpp_vec - ocpn_vec;

figure(4)
yyaxis left
plot(soc, ocv_vec, 'k-'); hold on
plot(soc, ocpp_vec, 'b-')
ylabel('Voltage[V]')
yyaxis right
plot(soc, ocpn_vec, 'r')
xline(1,'k--')
xline(0,'k--')
xlim([-0.2,1.3]);
ylabel('Voltage')
title('OCV from OCP')
xlabel('SOC')
 set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')

 % save('OCP2OCV')