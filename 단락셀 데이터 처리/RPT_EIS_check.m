
clear; clc; close all;

folder = 'G:\공유 드라이브\BSL-Data\카이스트_단락셀\카이스트 단락셀\2차 셀 데이터\RPT\RPT(edit)';
file = 'RPT_test_2_02_EIS.txt';

% 데이터 불러오기
opts = detectImportOptions([folder filesep file], 'NumHeaderLines', 13, 'VariableNamingRule', 'preserve');

data = readtable([folder filesep file], opts);

% 필요한 열 추출

time_cell = data.Var2;    % 두 번째 열
z_re = data.Var12; % 열 두 번째 열
z_im = data.Var13;% 열 세 번째 열
cycle = data.Var5;
color = lines(max(cycle));

for i = 1:max(cycle);

% 원하는 사이클만 골라 사용

idx = (cycle == i);

figure(i)
plot(z_re(idx),-z_im(idx),'o','color',color(i,:),'DisplayName', ['EIS' num2str(i)] );
xlabel("Z'")
ylabel("z''")
title ([file])
legend('show');
grid on;
box on;

end