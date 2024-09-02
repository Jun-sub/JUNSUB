% 충전 용량, 방전 용량 동시
clear; clc; close all;

folder = 'G:\공유 드라이브\BSL-Data\카이스트_단락셀\카이스트 단락셀\2차 셀 데이터\Aging\사이클 데이터';
files = dir(folder);

tic; 
% 모든 파일에 대해 반복
for k = 3:length(files)
    FileName = files(k).name;
    data_path = fullfile(folder, FileName);
    
    data = readtable(data_path, "FileType","text",'NumHeaderLines', 11, 'ReadVariableNames', false);
    
    % 두 번째 컬럼이 사이클, 세 번째 컬럼이 충전 용량, 네 번째 컬럼이 방전 용량
    cycle = data.Var2;
    charge_capacity = data.Var3;
    discharge_capacity = data.Var4;

    % 사이클 번호
    unique_cycles = unique(cycle);
    
    % Figure 초기화
    figure;
    hold on;
    
    % 사이클별로 그래프 그리기
    skip_flag = false; % 이상 데이터가 발견되면 다음 사이클로 넘어갈지 결정
    for i = 1:length(unique_cycles)-2
        current_cycle = unique_cycles(i);
        
        %이상 충전 데이터 제외
        % if i < length(unique_cycles) % 마지막 사이클에는 다음 사이클이 없으므로 범위 체크
        %     next_cycle = unique_cycles(i + 1);
        %     if charge_capacity(cycle == current_cycle(end))*1.05 < charge_capacity(cycle == next_cycle(end))
        %         fprintf('%d cycle abnormal charging data detected at cycle %d. Skipping subsequent cycles.\n', k, current_cycle);
        %         skip_flag = true;
        %         break; % 현재 파일의 이후 사이클을 모두 건너뛰기 위해 break
        %     end
        % end
        % 
        % 현재 사이클에 해당하는 인덱스 찾기
        idx = cycle == current_cycle;
       
        % 그래프 그리기
        plot(cycle(idx), charge_capacity(idx), 'bo', 'MarkerFaceColor','b');
        plot(cycle(idx), discharge_capacity(idx), 'ro', 'MarkerFaceColor','r');
    end

    % 사이클 건너뛸 경우에 대한 처리
    if skip_flag
        fprintf('Skipping remaining cycles in file %s.\n', FileName);
    end
    
    xlabel('Cycle Number');
    ylabel('Capacity [Ah]');
    title(strjoin(strsplit(FileName,'_'),' '));
    grid on;
    box on;
    xlim([min(cycle)-1, max(cycle)])
    xticks(min(cycle)-1:1:max(cycle));
    hold off;
    legend('Charge','Discharge','Location','Southwest');
   
end
toc;