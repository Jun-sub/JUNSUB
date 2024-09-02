clear; clc; close all;

folder = 'G:\공유 드라이브\BSL-Data\카이스트_단락셀\카이스트 단락셀\2차 셀 데이터\RPT\RPT(edit)';
file = ['240710_RPT_protocol_test_1_01_DC.txt'];

% 데이터 불러오기

data = readtable([folder filesep file], 'FileType', 'text', 'NumHeaderLines', 24, 'ReadVariableNames',false);

% 필요한 열 추출
time_cell = data.Var2;    % 두 번째 열
current = data.Var9; % 아홉 번째 열
voltage = data.Var10;% 열 번째 열

% cell 배열을 duration 형식으로 변환
cycle = data.Var5; % 다섯 번째 열이 사이클 정보
time = duration(time_cell,'InputFormat','hh:mm');

for i = 1:max(cycle);

% 원하는 사이클만 골라 사용

idx = (cycle == i);

figure()
hold on

yyaxis left
plot(time(idx), current(idx), 'b-', 'LineWidth', 1, 'DisplayName', 'Current(A)');
ylabel('Current (A)');

yyaxis right
plot(time(idx), voltage(idx), 'r-', 'LineWidth', 1, 'DisplayName', 'Voltage(V)');
ylabel('Voltage (V)');

xlabel('Time (s)');
title(['Time vs Current/Voltage for ', file]);
legend(['cycle', num2str(i)]);
grid on;
box on;

end

figure()
yyaxis left
plot(time, current, 'b-', 'LineWidth', 1, 'DisplayName', 'Current(A)');
ylabel('Current (A)');

yyaxis right
plot(time, voltage, 'r-', 'LineWidth', 1, 'DisplayName', 'Voltage(V)');
ylabel('Voltage (V)');

xlabel('Time (s)');
title(['Time vs Current/Voltage for ', file]);
legend('show');
grid on;
box on;