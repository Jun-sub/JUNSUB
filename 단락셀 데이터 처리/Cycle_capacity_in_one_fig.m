clc; clear; close all;

folder_path = 'G:\공유 드라이브\BSL-Data\카이스트_단락셀\카이스트 단락셀\3차 셀 데이터\Ref_txt';

folder = dir(folder_path);


cmap = jet(20); %colormap for dot
 cline = [0.8 0 0.2];


figure
hold on;
tic
for i = 3:5 %length(folder)
    data_path = fullfile(folder_path, folder(i).name);
    data = readtable(data_path,"FileType","text", 'NumHeaderLines',11,'ReadVariableNames',false);
    
    % if mod(i-2,6) == 1 %다섯 사이클마다 색 한번씩 변환
    %         cline = [rand(1) 0 rand(1)*0.3]; %매 사이클의 색상 정도 RGB
    %         fprintf('Count %d \n', floor((i-2)/6))
    % end
    
    cycle = data.Var2;
    cycle = cycle(1:end-2); %마지막과 그 전 데이터는 제외
    discharge_capacity = data.Var4;

    plot(cycle, discharge_capacity(cycle), '-', 'LineWidth', 1.5,'Color',cline);
    plot_name{i-2} = strjoin(strsplit(folder(i).name(1:end-4),'_'),' '); % 범례 지정용


    %마지막 사이클에만 마커 추가
    cycle_fin(i-2) = max(cycle);
    capa_fin(i-2) = discharge_capacity(max(cycle));
    
end
hold off;

title('Cycle capacity','FontWeight','bold')
xlabel('Cycle','FontWeight','bold')
ylabel('Calacity[Ah]','FontWeight','bold')
grid on;
box on;
xlim([1 20])
set(gca,'Fontunits','points','FontSize',20,'FontName','Times New Roman','PlotBoxAspectRatio',[1 1 1])
% legend(plot_name)

figure(1)
hold on;
for k = 1:length(cycle_fin)
    scatter(cycle_fin(k),capa_fin(k),50,'bo','filled','HandleVisibility','off');
end
hold off

% cbar = colorbar;
% colormap(jet(20));
% clim([0,20]);
% cbar.Label.String = 'Cycle';
% cbar.Label.FontWeight = 'bold';
% cbar.Label.FontSize = 10;

toc
