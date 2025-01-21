% This code includes
% 1) OCP, OCV initial plot 
% 2) OCV model fitting to obtain stoichiometry (para)
% 3) dV/dQ calculation and plot

clear; clc; close all


%% Configuration

% folder = 'G:\공유 드라이브\BSL-Data\LGES\2차 실험\OCP\데이터 변환';
% filename_ocpn = 'Processed_Data_AHC.mat';
% filename_ocpp = 'Processed_Data_CHC.mat';
% filename_ocv = 'Processed_Data_OCV.mat';


folder = 'G:\공유 드라이브\EE6211_2024\Data\pOCV_and_pOCP';

filename_ocpn = 'AHC_(5)_OCV_C20.mat';
filename_ocpp = 'CHC_(5)_OCV_C20.mat';
filename_ocv = 'FCC_(5)_OCV_C20.mat';
%smooth degree & OCV direction set
num_move = 100;
direc = 1; % 0: charge, 1: discharge

%stoi initial
x0 = 0.0255;
x1 = 0.8900;
y0 = 0.8568;
y1 = 0.2128;

para_ini = [x0 x1 y0 y1];
%soc
soc_01 = 0:0.001:1;
x_soc01 = x0 + (x1-x0)*soc_01;
y_soc01 = y0 + (y1-y0)*soc_01;

%data load set

if direc == 0;
    p_direc = 'chg';
    n_direc = 'dis';
    o_direc = 'chg';
   
elseif direc == 1;
    p_direc = 'dis';
    n_direc = 'chg';
    o_direc = 'dis';
else 
    error('Assign 0:charge or 1: discharge in direc')
end

%data load
ocpn = load([folder filesep filename_ocpn]);
ocpp = load([folder filesep filename_ocpp]);

ocpn = ocpn.OCV_golden.(['OCV' n_direc]); %ocp
ocpp = ocpp.OCV_golden.(['OCV' p_direc]);

ocv = load([folder filesep filename_ocv]);
ocv = ocv.OCV_golden.(['OCV' p_direc]);

%% OCP, OCV
%OCV plot
figure()
t = tiledlayout(3,1,"TileSpacing","compact",'Padding','compact');

nexttile
plot(ocv(:,1),ocv(:,2),'k-');
title('OCV')

nexttile
plot(ocpn(:,1),ocpn(:,2),'b-');
title('OCPn');

nexttile
plot(ocpp(:,1),ocpp(:,2),'r-');
title('OCPp')


%OCV Calc
ocpn_calc = interp1(ocpn(:,1),ocpn(:,2),x_soc01);
ocpp_calc = interp1(ocpp(:,1),ocpp(:,2),y_soc01);
ocv_calc = ocpp_calc - ocpn_calc;

% calc plot check
figure()
hold on;
%ocpn
yyaxis('left')
ylabel('Voltage [V]')
plot(soc_01,ocpn_calc,'b-'); 

%ocv & ocpp
yyaxis("right")
plot(ocv(:,1),ocv(:,2),'k-'); 
plot(soc_01,ocv_calc,'g-');
plot(soc_01,ocpp_calc,'r-');
ylabel('Voltage [V]')

hold off;

set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman','XGrid','on','YGrid','on','XLim',[min(soc_01) max(soc_01)]);
legend('OCPn','OCV','OCV calc','OCPp')

xlabel('soc with stoi')

%% OCV fitting

lb = [0 0.9 0.8 0.15];
% lb(2) = 0.8;
ub = [0.1 1 0.95 0.25];

options = optimoptions("fmincon","Display","iter","MaxIterations",100);

fhandle_cost = @(para)func_cost(ocv(:,1),para,ocpp,ocpn,ocv(:,2));

para_hat = fmincon(fhandle_cost,para_ini,[],[],[],[],lb,ub,[],options);

% opimized OCV plot
ocv_hat = func_ocv_model(soc_01,para_hat,ocpp,ocpn);

figure()
hold on;
plot(soc_01, ocv_calc); hold on
plot(soc_01, ocv_hat); hold on
plot(ocv(1:50:end,1),ocv(1:50:end,2),'o')
legend({'initial gues','model fit','data'})
pause(0.1)
hold off


% OCV with extended stoi
soc_ex = -0.2:0.01:1.2;
x_soc_ex = para_hat(1) + (para_hat(2)-para_hat(1))*soc_ex;
y_soc_ex = para_hat(3) + (para_hat(4)-para_hat(3))*soc_ex;

%OCV Calc
ocpn_calc_ex = interp1(ocpn(:,1),ocpn(:,2),x_soc_ex);
ocpp_calc_ex = interp1(ocpp(:,1),ocpp(:,2),y_soc_ex);
ocv_calc_ex = ocpp_calc_ex - ocpn_calc_ex;

% calc plot check
figure()
hold on;
%ocpn
yyaxis('left')
ylabel('Voltage [V]')
plot(soc_ex,ocpn_calc_ex,'b-'); 

%ocv & ocpp
yyaxis("right")
plot(soc_01,ocv_hat,'g-');
plot(soc_ex,ocpp_calc_ex,'r-');
ylabel('Voltage [V]')

hold off;
set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman','XGrid','on','YGrid','on','XLim',[min(soc_ex) max(soc_ex)]);

xlabel('soc')

%% dV/dQ Calculation

%dVdQ calc
if direc == 0
    dVdQ_ocv = diff(ocv(:,2))./diff(ocv(:,3));
    dVdQ_ocpp = diff(ocpp(:,2))./diff(ocpp(:,3));
    dVdQ_ocpn = diff(ocpn(:,2))./diff(ocpn(:,3));
elseif direc == 1
     dVdQ_ocv = -diff(ocv(:,2))./diff(ocv(:,3));
    dVdQ_ocpp = -diff(ocpp(:,2))./diff(ocpp(:,3));
    dVdQ_ocpn = -diff(ocpn(:,2))./diff(ocpn(:,3));
end

%Smooth
dVdQ_ocv = movmean(dVdQ_ocv,num_move);
dVdQ_ocpp = movmean(dVdQ_ocpp,num_move);
dVdQ_ocpn = movmean(dVdQ_ocpn,num_move);

%Visuallization

x_soc = linspace(-0.2,1.2,length(dVdQ_ocpn));
y_soc = linspace(-0.2,1.2,length(dVdQ_ocpp));

x_stoi = para_hat(1) + (para_hat(2)-para_hat(1))*x_soc;
y_stoi = para_hat(3) + (para_hat(4)-para_hat(3))*y_soc;

figure
hold on;
plot(ocv(1:end-1,1),dVdQ_ocv,'k-','LineWidth',1);
% plot(1-ocpp(1:end-1,1),dVdQ_ocpp,'r:','LineWidth',2);
% plot(ocpn(1:end-1,1),dVdQ_ocpn,'b:','LineWidth',2);

plot(y_stoi,dVdQ_ocpp,'r:','LineWidth',1);
plot(x_stoi,dVdQ_ocpn,'b:','LineWidth',1);

grid on;
box on;
xlabel('SOC')
ylabel('dVdQ')
legend('OCV','OCPp','OCPn')
pbaspect([1 1 2]); % X:Y 비율을 1:1로 설정 (필요에 따라 조정)
ylim_top = 2*max(dVdQ_ocv((ocv(:,1) > 0.2) & (ocv(:,1) < 0.8)));
ylim([0 ylim_top])

%% Function

function cost = func_cost(soc,para,ocpp,ocpn,ocv_data)

ocv_model = func_ocv_model(soc,para,ocpp,ocpn);

cost = sum(sqrt((ocv_data - ocv_model).^2));

end

function ocv_model = func_ocv_model(soc,para,ocpp,ocpn)

x0 = para(1);
x1 = para(2);
y0 = para(3);
y1 = para(4);

x_soc01 = x0 + (x1-x0)*soc;
y_soc01 = y0 + (y1-y0)*soc;

ocpp_soc01 = interp1(ocpp(:,1),ocpp(:,2),y_soc01,'linear','extrap');
ocpn_soc01 = interp1(ocpn(:,1),ocpn(:,2),x_soc01,'linear','extrap');

ocv_model = ocpp_soc01 - ocpn_soc01;

end