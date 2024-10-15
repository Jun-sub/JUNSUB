clc;clear;close all

%% VERSION SUMMARY
% 24-10-15
% 1) Use ocpn_chg 
% 2) Moving Average: i. moving average in ocv, not in dvdq (also have speed-up effect)
% 3) Moving Average: ii. define by number of points 
% 4) Weighting: i. relative rmse is used (diverging part is naturally de-weighted). Deleted codes finding the minimum dvdq points.

%% Configuration

n_points = 95; %for soc range resolution
direc = 0; % 0 for chg direction, 1 for discharge direction



% charge or discharge assign
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

% load OCP datas
folder = 'G:\공유 드라이브\BSL-Data\LGES\2차 실험\OCP\데이터 변환';

filename_ocpn = 'Processed_Data_AHC.mat';
filename_ocpp = 'Processed_Data_CHC.mat';
filename_ocv = 'Processed_Data_OCV.mat';

%ocpn
ocpn = load([folder filesep filename_ocpn]);
ocpn_raw = ocpn.OCV_golden.(['OCV' n_direc]); %ocpn

x_raw = ocpn_raw(:,1); %soc of ocpn
x = linspace(min(x_raw),max(x_raw),n_points)'; %soc n_points만큼 나눠서 표현
ocpn_raw = ocpn_raw(:,2); %ocpn
ocpn_mva = movmean(ocpn_raw,round(length(ocpn_raw)/n_points)); %smoothed ocpn
ocpn = interp1(x_raw,ocpn_mva,x);

ocpn = [x ocpn]; %movmean ocpn inter1 한 데이터

clear x_raw x ocpn_raw ocpn_mva

%ocpp
ocpp = load([folder filesep filename_ocpp]);
ocpp_raw = ocpp.OCV_golden.(['OCV' p_direc]);

y_raw = ocpp_raw(:,1); %soc of ocpp
y = linspace(min(y_raw),max(y_raw),n_points)'; %soc n_points만큼 나눠서 표현
ocpp_raw = ocpp_raw(:,2); %ocpn
ocpp_mva = movmean(ocpp_raw,round(length(ocpp_raw)/n_points)); %smoothed ocpn
ocpp = interp1(y_raw,ocpp_mva,y);

ocpp = [y ocpp]; %movmean ocpp inter1 한 데이터

clear y_raw y ocpp_raw ocpp_mva

%ocv
ocv = load([folder filesep filename_ocv]);
ocv_raw = ocv.OCV_golden.(['OCV' p_direc]);

