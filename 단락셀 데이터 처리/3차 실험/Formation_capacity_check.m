clc; clear; close all;

folder_path = 'G:\공유 드라이브\BSL-Data\카이스트_단락셀\카이스트 단락셀\3차 셀 데이터\본 실험\Formation';

folder = dir(folder_path);
folder = folder(3:end,:);

% sort name sequence
name = {folder.name};

numbers = cellfun(@(x) sscanf(x,'%d'),name);
[sorted_num,idx] = sort(numbers);


figure
hold on;
for i = 1:length(folder)
    
    file_name{i} = folder(idx(i)).name;

    data_now = readmatrix(fullfile(folder_path,file_name{i}),"NumHeaderLines",9);

    capacity(i,1) = data_now(end,4);
    
end

scatter(1:length(folder),capacity,40,'o','filled','MarkerEdgeColor','black');

grid on;
box on;
xlabel('Sample number')
ylabel('Ah')
set(gca,"PlotBoxAspectRatio",[1 1 1])
hold off;