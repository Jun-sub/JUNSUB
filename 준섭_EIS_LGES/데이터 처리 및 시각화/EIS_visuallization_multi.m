clear; clc; close all;

folder = 'G:\공유 드라이브\BSL-Data\LGES\2차 실험\RPT';
file = 'CHC_C_20_RPT 2_02_EIS.txt';

% 데이터 불러오기
opts = detectImportOptions([folder filesep file], 'NumHeaderLines', 13, 'VariableNamingRule', 'preserve');
data = readtable([folder filesep file], opts);

% 필요한 열 추출
time_cell = data.Var2;    % 두 번째 열
z_re = data.Var12;        % 열 두 번째 열
z_im = data.Var13;        % 열 세 번째 열
cycle = data.Var5;

% 사이클 최대값 구하기
max_cycle = max(cycle);

% ColorOrder 설정: 다양한 색상을 사용하기 위해 미리 정의
color_order = lines(max_cycle);  % 예시로 MATLAB 기본 색상 사용
figure
    hold on
% 각 사이클마다 플롯
for i = 2:max_cycle
    idx = cycle == i;

    
    plot(z_re(idx), -z_im(idx), 'o', 'Color', color_order(i, :), 'DisplayName', ['Cycle ', num2str(i)]);
    xlabel("Z'")
    ylabel("Z''")
    title(file)
    legend('show')
    grid on
    box on
    
end

hold off