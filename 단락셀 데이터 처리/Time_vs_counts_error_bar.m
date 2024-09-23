clear; clc; close all;

%x axis: DOE, y axis: time [s]
folder = 'G:\공유 드라이브\BSL-Data\카이스트_단락셀\카이스트 단락셀\2차 셀 데이터\Aging';
files = dir(folder);

files = files(3:end-1);

data_st = struct('name',[],'time',[],'index',[]);
n = 1;
for i = 1:length(files)
% clear time_sum; clear Q_sum
    tic_cycle = tic;
    
    file_name = files(i).name;
    data_now = readtable([folder filesep file_name], "FileType","text",'NumHeaderLines', 20, 'ReadVariableNames', false);
    
    % data_trim
    data_now = data_now(data_now.Var18 ~= 0,:);

    time = data_now.Var2;
    cycle = data_now.Var3;
    
    name_idx = strfind(file_name,'_');
    name = file_name(1:name_idx(2) - 1);

    time_end = time(cycle == max(cycle) - 2);
    time_end = seconds(time_end(end));
   
    %data_st assign
    data_st(i).name = name;
    data_st(i).time = time_end;
    
    %index assign
    

    fprintf('%d step of %d is processed \n',i,length(files))
    toc(tic_cycle);
end
  
for k = 1:length(files)-1 
data_st(k).index = n;
    doe_f = split(data_st(k).name,'_');
    doe_s = split(data_st(k+1).name,'_');

    if ~strcmp(doe_f{1}, doe_s{1})
       n = n+1;
    end

    data_st(k+1).index = n;
end

%% Calculate mean and mark error bar for each DOE
data_sum = [];
for j = 1:data_st(end).index
   data_cal = data_st([data_st.index] == j); %Bring specific index
    for f = 1:length([data_cal.index])
        time_doe(f) = data_cal(f).time;
        name = split(data_cal(f).name,'_');
        name = name{1};
    end
    mean_data = mean(time_doe);
    std_data = std(time_doe);
    
    data_sum(j).name = name;
    data_sum(j).mean = mean_data;
    data_sum(j).std = std_data;
end

%plot
figure(1)
errorbar([data_sum.mean], [data_sum.std], 'ro', 'MarkerSize',6,'LineWidth',1.5, 'MarkerFaceColor','r')

%x-axis set
set(gca,'XTick',1:length({data_sum.name}));
set(gca,'XTickLabel',{data_sum.name});
xlim([0 length({data_sum.name})+1]);


xtickangle(45)
xlabel('DOE');
ylabel('Max time [s]');

grid on;
box on;
