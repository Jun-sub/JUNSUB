clear; clc; close all


%% Configuration
%data load

folder = 'G:\공유 드라이브\BSL-Data\LGES\2차 실험\OCP\데이터 변환';
filename_ocpn = 'Processed_Data_AHC.mat';
filename_ocpp = 'Processed_Data_CHC.mat';
filename_ocv = 'Processed_Data_OCV.mat';

ocpn = load([folder filesep filename_ocpn]);
ocpp = load([folder filesep filename_ocpp]);

ocpn = ocpn.OCV_golden.OCVdis;
ocpp = ocpp.OCV_golden.OCVchg;

ocv = load([folder filesep filename_ocv]);
ocv = ocv.OCV_golden.OCVchg;

%% dV/dQ Calculation
%data load

% soc_ocv = linspace(min(ocv(:,1)),max(ocv(:,1)),round(length(ocv(:,1)))/3);
% V_ocv = ocv(:,2);
% Q_ocv = ocv(:,3);
% 
% % V_model = ocv_hat
% % Q_model =
% 
% soc_ocpp = linspace(min(ocpp(:,1)),max(ocpp(:,1)),round(length(ocpp(:,1)))/3);
% V_ocpp = ocpp(:,2);
% Q_ocpp = ocpp(:,3);
% 
% soc_ocpn = linspace(min(ocpn(:,1)),max(ocpn(:,1)),round(length(ocpn(:,1)))/3);
% V_ocpn = ocpn(:,2);
% Q_ocpn = ocpn(:,3);

%dVdQ calc
dVdQ_ocv = diff(ocv(:,2))./diff(ocv(:,3));
dVdQ_ocv = movmean(dVdQ_ocv,100);

dVdQ_ocpp = diff(ocpp(:,2))./diff(ocpp(:,3));
dVdQ_ocpp = movmean(dVdQ_ocpp,100);

dVdQ_ocpn = diff(ocpn(:,2))./diff(ocpn(:,3));
dVdQ_ocpn = movmean(dVdQ_ocpn,100);
%Visuallization

figure
hold on;
plot(ocv(1:end-1,1),dVdQ_ocv,'k-','LineWidth',2);
plot(ocpp(1:end-1,1),dVdQ_ocpp,'r:','LineWidth',2);
plot(ocpn(1:end-1,1),dVdQ_ocpn,'b:','LineWidth',2);

grid on;
box on;
xlabel('SOC')
ylabel('dVdQ')
legend('OCV','OCPp','OCPn')
pbaspect([1 1 2]); % X:Y 비율을 1:1로 설정 (필요에 따라 조정)
ylim_top = 2*max(dVdQ_ocv((ocv(:,1) > 0.2) & (ocv(:,1) < 0.8)));
ylim([0 ylim_top])