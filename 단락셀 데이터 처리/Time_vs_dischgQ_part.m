clear; clc; close all;

%순서: 0인 데이터 trim --> 데이터 스텝별 나누기 --> 시간 플랏
folder = 'G:\공유 드라이브\BSL-Data\카이스트_단락셀\카이스트 단락셀\2차 셀 데이터\Aging';
files = dir(folder);

files = files(14:18,:);  %부분별 플랏할 때 사용
c_map = jet(length(files));
cline = [0.8 0 0.2];

tic_overall = tic;
figure() %time vs Q
hold on;
for i = 1:length(files)
clear time_sum; clear Q_sum
    tic_cycle = tic;
    
    file_name = files(i).name;
    data = readtable([folder filesep file_name], "FileType","text",'NumHeaderLines', 20, 'ReadVariableNames', false);
    
    % data_trim
    data = data(data.Var18 ~= 0,:);

    time = data.Var2;
    cycle = data.Var3;
    Q = data.Var18;
    
    % 마지막 강제사이클 제외
    for k = 1:max(cycle)-2
        cycle_idx = (cycle == k);
        time_end = time(cycle_idx);
        Q_end = Q(cycle_idx);
        
        time_sum(k) = time_end(end);
        Q_sum(k) = Q_end(end);
 
    end
    plot(time_sum, Q_sum, 'k-','LineStyle','-','Color',cline,'LineWidth',1)
    scatter(time_sum(end),Q_sum(end),30,'bo','filled')

    legend_name(i) = {file_name};

    data_t_pdf(i) = time_sum(end);
    data_c_pdf(i) = max(cycle)-2;
    fprintf('%d step of %d is processed \n',i,length(files))
    toc(tic_cycle);
end
hold off;

title('Time vs Q');
xlabel('Time [s]');
ylabel('Capacity [mAh]');
legend(legend_name(i));
grid on; box on;
set(gca,'Fontunits','points','FontSize',10,'FontName','Times New Roman','PlotBoxAspectRatio',[1 1 1])

% Normalized time dist
data_t_pdf = data_t_pdf(1:end);
data_t_pdf = seconds(data_t_pdf);

figure() %time vs probability density bar
histogram(data_t_pdf,'NumBins', 10,'Normalization','count');
x = linspace(min(data_t_pdf), max(data_t_pdf), 100);
pdf = normpdf(x,mean(data_t_pdf), std(data_t_pdf));

hold on;
plot(x, pdf * numel(data_t_pdf) * (max(data_t_pdf) - min(data_t_pdf)) / 10, 'r', 'LineWidth', 2); % 빨간색 선으로 PDF 표시, 이부분 이해 안됨.
xlabel('Time [s]');
ylabel('counts');
title('Time NDPDF');
legend('Histogram', 'ND-PDF');
hold off;

% Cycle dist
data_c_pdf = data_c_pdf(1:end);

figure() %cycle vs probability density bar
histogram(data_c_pdf,'NumBins', 10,'Normalization','count','FaceColor','y');
x = linspace(min(data_c_pdf), max(data_c_pdf), 100);
pdf = normpdf(x,mean(data_c_pdf), std(data_c_pdf));

hold on;
plot(x, pdf * numel(data_c_pdf) * (max(data_c_pdf) - min(data_c_pdf)) / 10, 'r', 'LineWidth', 2); % 빨간색 선으로 PDF 표시
xlabel('Cycle [n]');
ylabel('counts');
title('Cycle NDPDF');
legend('Histogram', 'ND-PDF');
hold off;
fprintf('전체')
toc(tic_overall);