clc;clear;close all

%% OCV fitting to OCP data to acquire stoichiometry infomation
%% [x0 x1 y0 y1]

%% Configuration

n_points = 400; %for soc range resolution
direc = 1; % 0 for chg direction, 1 for discharge direction



% charge or discharge assign
if direc == 0
    p_direc = 'chg';
    n_direc = 'dis';
    o_direc = 'chg';
   
elseif direc == 1
    p_direc = 'dis';
    n_direc = 'chg';
    o_direc = 'dis';
else 
    error('Assign 0:charge or 1: discharge in direc')
end

%% load OCP datas
folder = 'G:\공유 드라이브\BSL-Data\LGES\2차 실험\OCP\데이터 변환';

% folder = 'G:\공유 드라이브\EE6211_2024\Data\pOCV_and_pOCP';

filename_ocpn = 'Processed_Data_AHC.mat';
filename_ocpp = 'Processed_Data_CHC.mat';
filename_ocv = 'Processed_Data_OCV.mat';

% filename_ocpn = 'AHC_(5)_OCV_C20.mat';
% filename_ocpp = 'CHC_(5)_OCV_C20.mat';
% filename_ocv = 'FCC_(5)_OCV_C20.mat';

%ocpn
ocpn = load([folder filesep filename_ocpn]);
ocpn_raw = ocpn.OCV_golden.(['OCV' n_direc]); %ocpn

x_raw = ocpn_raw(:,1); %soc of ocpn
x = linspace(min(x_raw),max(x_raw),n_points)'; %soc n_points만큼 나눠서 표현
Qx_raw = ocpn_raw(:,3); %Q
ocpn_raw = ocpn_raw(:,2); %ocpn
ocpn_mva = movmean(ocpn_raw,round(length(ocpn_raw)/n_points)); %smoothed ocpn
ocpn = interp1(x_raw, ocpn_mva,x);

Qx = interp1(x_raw, Qx_raw, x, 'linear','extrap');

ocpn = [x ocpn Qx]; %anode stoi, movmean ocpn inter1, anode capacity


clear x_raw x ocpn_raw ocpn_mva

%ocpp
ocpp = load([folder filesep filename_ocpp]);
ocpp_raw = ocpp.OCV_golden.(['OCV' p_direc]);

y_raw = ocpp_raw(:,1); %soc of ocpp
y = linspace(min(y_raw),max(y_raw),n_points)'; %soc n_points만큼 나눠서 표현
Qy_raw = ocpp_raw(:,3); %Q
ocpp_raw = ocpp_raw(:,2); %ocpn
ocpp_mva = movmean(ocpp_raw,round(length(ocpp_raw)/n_points)); %smoothed ocpn
ocpp = interp1(y_raw,ocpp_mva,y);

Qy = interp1(y_raw, Qy_raw, y, 'linear','extrap');

ocpp = [y ocpp Qy]; %cathode stoi, movmean ocpp inter1, cathode capacity


clear y_raw y ocpp_raw ocpp_mva

%ocv
ocv = load([folder filesep filename_ocv]);
ocv_raw = ocv.OCV_golden.(['OCV' p_direc]);

z_raw = ocv_raw(:,1); %soc of ocpp
soc_01 = linspace(min(z_raw),max(z_raw),n_points)'; %soc n_points만큼 나눠서 표현
Qz_raw = ocv_raw(:,3); %Q
ocv_raw = ocv_raw(:,2); %ocpn
ocv_mva = movmean(ocv_raw,round(length(ocv_raw)/n_points)); %smoothed ocpn
ocv = interp1(z_raw,ocv_mva,soc_01);

Qz = interp1(z_raw, Qz_raw, soc_01, 'linear','extrap');

ocv = [soc_01 ocv Qz]; % z = soc

% plot ocpn, ocpp, ocv data
figure

t = tiledlayout(3,1,"TileSpacing","compact","Padding","compact");

nexttile
plot(ocpn(:,1),ocpn(:,2),'r-')
legend('ocpn raw')

nexttile
plot(ocpp(:,1),ocpp(:,2),'b-')
legend('ocpp raw')

nexttile
plot(ocv(:,1),ocv(:,2),'k-')
legend('ocv raw')

%% Fitting OCV curve
% ocv data generation from ocp
para_ini = [0.02, 0.925, 0.8067, 0.2360]; %x0 x1 y0 y1


% ocv_model0 = func_ocv_model(para_ini, soc_01, ocpn(:,1:2),ocpp(:,1:2));

% fitting to ocv exp data
function_handle = @(para,soc) func_ocv_model(para,soc,ocpn,ocpp);

options = optimoptions("lsqcurvefit",'Display','iter','MaxIterations',100,'FunctionTolerance',1e-9,"ConstraintTolerance",1e-9);

[para_hat, resnorm] = lsqcurvefit(function_handle, para_ini, ocv(:,1), ocv(:,2),[],[],options); %[0 0.9 0.75 0.15],[0.1 1 0.9 0.25],options)

% fitting data plot
ocv_model1 = func_ocv_model(para_hat, soc_01, ocpn(:,1:2),ocpp(:,1:2));

figure
plot(ocv(:,1),ocv(:,2),'k-',soc_01(1:4:end),ocv_model1(1:4:end),'bo','MarkerSize',5)
legend('Exp data','Model OCV','Location','best')

% Stoi range extension
soc_ex = linspace(-0.2,1.2,1.5*n_points);

x0 = para_hat(1);
x1 = para_hat(2);
y0 = para_hat(3);
y1 = para_hat(4);

soc_hat_x = x0 + (x1-x0)*soc_01;
soc_hat_y = y0 + (y1-y0)*soc_01;

soc_ex_x = x0 + (x1-x0)*soc_ex;
soc_ex_y = y0 + (y1-y0)*soc_ex;

ocpn_ex = interp1(ocpn(:,1),ocpn(:,2),soc_ex_x);
ocpp_ex = interp1(ocpp(:,1),ocpp(:,2),soc_ex_y);

figure
hold on;
yyaxis left
plot(soc_01,ocv_model1,'k-',soc_ex, ocpp_ex,'b-')

yyaxis right
plot(soc_ex, ocpn_ex,'r-')
legend('Exp data','ocpp','ocpn','Location','best')
hold off


%% dVdQ Calculation
%dV
dVn = diff(ocpn(:,2));
dVp = diff(ocpp(:,2));
dVo = diff(ocv_model1);

%dQ
dQn = diff(flip(Qz));
dQp = diff(flip(Qz));
% dQn = diff(Qy);
% dQp = diff(Qx);
dQo = diff(Qz);

%plot
if direc == 0
    figure
    hold on
    plot(soc_01(1:end-1),dVo./dQo,'k-'); %dVdQ ocv
    plot(soc_01(1:end-1),dVn./dQn,'r-'); %dVdQ ocpn
    plot(soc_01(1:end-1),flip(dVp./dQp),'b-'); %dVdQ ocpp
    hold off
legend('OCV','OCPn','OCPp')
elseif direc == 1
    figure
    hold on
    plot(soc_01(1:end-1),-dVo./dQo,'k-'); %dVdQ ocv
    plot(soc_01(1:end-1),-dVn./dQn,'r-'); %dVdQ ocpn
    plot(soc_01(1:end-1),-flip(dVp./dQp),'b-'); %dVdQ ocpn
    hold off
    legend('OCV','OCPn','OCPp')
end 

dvdqo = abs(dVo./dQo);
y_top = 2*(max(dvdqo((soc_01 > 0.2) & (soc_01 < 0.8))));
ylim([0 y_top])
function [ocv_model] = func_ocv_model(para,soc,ocpn,ocpp)
    
    x0 = para(1);
    x1 = para(2);
    y0 = para(3);
    y1 = para(4);

    x = x0 + (x1-x0)*soc;
    y = y0 + (y1-y0)*soc;

    calc_ocpn = interp1(ocpn(:,1),ocpn(:,2),x,"linear","extrap");
    calc_ocpp = interp1(ocpp(:,1),ocpp(:,2),y,"linear","extrap");

    ocv_model = calc_ocpp - calc_ocpn;
end
