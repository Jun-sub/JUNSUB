clear; clc; close all;

header_row = 13;

folder_path = 'G:\공유 드라이브\BSL-Data\카이스트_단락셀\카이스트 단락셀\2차 셀 데이터\RPT\RPT(edit)';
folder = dir(folder_path)

for i = 1:length(folder)
    
    file = folder(i).name;
    
    if ~contains(file,'EIS')
        continue;
    end



% 데이터 불러오기
data = readtable([folder filesep file], 'NumHeaderLines', header_row, 'ReadVariableNames',false);

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

end
