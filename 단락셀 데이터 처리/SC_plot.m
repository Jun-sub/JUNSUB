%단락셀 데이터를 불러와서 이어서 플랏
clc;, clear all;, close;

% Sample: 
% 폴더 지정하면 각 파일 분류할 수 있게.
full_folder_path = 'G:\공유 드라이브\BSL-Data\카이스트_단락셀\카이스트 단락셀\3차 셀 데이터\Ref_txt'
full_folder = dir(full_folder_path);

for i = 3:5 %18
folder_path = ([full_folder_path filesep full_folder(i).name]);
folder = dir(fullfile(folder_path, '*.txt'));

data_form = readtable([folder_path filesep folder(2).name],"FileType","text",'NumHeaderLines', 11, 'ReadVariableNames',0);
data_ocv = readtable([folder_path filesep folder(3).name], "FileType","text",'NumHeaderLines', 11, 'ReadVariableNames',0);
data_age = readtable([folder_path filesep folder(1).name], "FileType","text",'NumHeaderLines', 11, 'ReadVariableNames',0);

data_RPT = readtable([folder_path filesep folder(4).name], "FileType", "text" , 'ReadVariableNames', 0);

% 데이터 시간열 추출 (Formation --> OCV --> EIS_DCIR --> Aging 순으로 진행)

t_form = seconds(data_form.Var2);
t_ocv = t_form(end,1) + seconds(data_ocv.Var2); % 이전 포메이션 시간 끝 + OCV 시간

% RPT 시간 데이터가 119 열 까지와 120 이후부터 형식이 다름.
t_RPT = t_ocv(end,1) + seconds(duration(data_RPT.Var2)); % OCV 끝나는 시간 + RPT 시간

t_age = t_RPT(end,1) + seconds(data_age.Var2); % RPT 끝나는 시간 + Aging 시간

% 데이터 전압 추출
v_form =data_form.Var8;
v_ocv = data_ocv.Var8;
v_age = data_age.Var8;
v_RPT = data_RPT.Var10;

% 데이터 전류 추출
i_form =data_form.Var7;
i_ocv = data_ocv.Var7;
i_age = data_age.Var7;

i_RPT = data_RPT.Var9;

% 그래프 플랏
figure(i)

yyaxis left
plot (t_form, v_form, 'r-'); hold on;
plot (t_ocv, v_ocv, 'r-'); hold on;
plot (t_RPT, v_RPT, 'r-'); hold on;
plot (t_age, v_age,'r-'); hold on;

set (gca, 'YColor', 'r');
ylim auto;
xlabel ('time [s]')
ylabel ('Voltage [V]')

yyaxis right
plot (t_form, i_form, 'b-'); hold on;
plot (t_ocv, i_ocv,'b-'); hold on;
plot (t_RPT, i_RPT,'b-'); hold on;
plot (t_age, i_age,'b-'); hold on;

set (gca, 'YColor', 'b');
xlabel ('time [s]')
ylabel ('Current [A]')
ylim auto
title(folder(1).name(14:18))
end