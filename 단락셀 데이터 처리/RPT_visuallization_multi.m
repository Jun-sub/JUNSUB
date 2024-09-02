%RPT visuallization

clear; clc; close all;

%데이터를 불러와서 이어서 플랏

% 폴더 지정하면 각 파일 분류할 수 있게.
full_folder_path = 'G:\공유 드라이브\BSL-Data\카이스트_단락셀\카이스트 단락셀\2차 셀 데이터\RPT\RPT(edit)'
full_folder = dir(full_folder_path);




for i = 3:length(full_folder)
    name = full_folder(i).name;
    path = fullfile(full_folder_path, name); 
   
    if ~contains(name, 'DC') || contains(name, 'EIS')
        continue; % DC.txt를 포함하지 않거나 EIS.txt를 포함하면 다음 반복으로 넘어감
    end

    %if max(ismember(strsplit(full_folder(i).name,'_'),'DC.txt'));

        data = readtable(path, 'NumHeaderLines', 23, ReadVariableNames=false);

    % time_str을 테이블에서 추출
    
    time = data.Var2; % 시간 문자열 열 추출
    I = data.Var9;
    V = data.Var10;
    cycle = data.Var5;

    
 
    %plot
    figure()
    hold on
    
    yyaxis left
    plot(time,I,'b-', 'LineWidth', 1)
    ylabel('Current [A]')
    
    
    yyaxis right
    plot(time,V,'r-', 'LineWidth', 1)
    ylabel('Voltage [V]')
    hold off
    
    set(gca,'Box','on',... %Axis Properties: BOX   
        'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
        'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    
    grid on;
    title (['I/V vs time graph',num2str(i)])
    xlabel('Time [s]')
    
    legend ('Current', 'Voltage', 'Location', 'southwest')
    
    %step 별 확대 플랏
    
    
    for k = 1:max(cycle);
    
    idx = (cycle == k);
    
    figure()
    hold on
    
    yyaxis left
    plot(time(idx), I(idx), 'b-', 'LineWidth', 1, 'DisplayName', 'Current(A)');
    ylabel('Current (A)');
    
    yyaxis right
    plot(time(idx), V(idx), 'r-', 'LineWidth', 1, 'DisplayName', 'Voltage(V)');
    ylabel('Voltage (V)');
    
    xlabel('Time (s)');
    title('Time vs Current/Voltage profile');
    legend(['step', num2str(k)], 'Location','southwest');
    grid on;
    box on;
    
    end 
    end
