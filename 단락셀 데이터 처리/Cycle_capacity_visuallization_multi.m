% 충전 용량, 방전 용량 동시
clear; clc; close all;

folder = 'G:\공유 드라이브\BSL-Data\카이스트_단락셀\카이스트 단락셀\3차 셀 데이터\Ref_cycle_txt';
files = dir(folder);

folder_2 = 'G:\공유 드라이브\BSL-Data\카이스트_단락셀\카이스트 단락셀\3차 셀 데이터\Ref_txt';
files_2 = dir(folder_2);

tic; 
% 모든 파일에 대해 반복
figure
t = tiledlayout(6,6,"TileSpacing","tight","Padding","tight");
for k = 3:length(files)
    nexttile
    FileName = files(k).name;
    FileName_2 = files_2(k).name;

    data_path = fullfile(folder, FileName);
    data_path_2 = fullfile(folder_2, FileName_2);
    
    data = readtable(data_path, "FileType","text",'NumHeaderLines', 11, 'ReadVariableNames', false);
    data_2 = readtable(data_path_2, "FileType","text",'NumHeaderLines', 11, 'ReadVariableNames', false);

    % 두 번째 컬럼이 사이클, 세 번째 컬럼이 충전 용량, 네 번째 컬럼이 방전 용량
    cycle = data.Var2;
    charge_capacity = data.Var3;
    discharge_capacity = data.Var4;
    last_cycle_time = data_2.Var2(end);
    
    % 사이클 번호
    unique_cycles = unique(cycle);
    
    % Figure 초기화
    % figure;
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
        % plot(cycle(idx), charge_capacity(idx), 'bo', 'MarkerFaceColor','b');
        plot(cycle(idx), discharge_capacity(idx), 'ro', 'MarkerFaceColor','r','MarkerSize',3);
        plot(seconds(last_cycle_time)/(24*60*60),discharge_capacity(end),'')
    end

    % 사이클 건너뛸 경우에 대한 처리
    if skip_flag
        fprintf('Skipping remaining cycles in file %s.\n', FileName);
    end
    
    xlabel('Cycle Number');
    ylabel('Capacity [Ah]');
    title_name = strsplit(FileName,'_');
    title(strjoin(title_name(1:2),'_'));
    grid on;
    box on;
    xlim([min(cycle)-2, max(cycle)+4])
    % xticks(min(cycle)-1:1:max(cycle));
    hold off;
    % legend('Charge','Discharge','Location','best');
    legend(['Discharge ',num2str(max(cycle))],[num2str(seconds(last_cycle_time)/(24*60*60)) ' day'],'Location','best');
    set(gca,'PlotBoxAspectRatio', [1 1 1])
end

set(gcf,"Position",[100 100 1200 1200])
toc;