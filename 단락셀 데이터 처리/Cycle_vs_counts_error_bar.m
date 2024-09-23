clear; clc; close all;

%순서: 0인 데이터 trim --> 데이터 DOE별 나누기(name,cycle,index) --> 각 DOE별 평균 및 error plot
folder = 'G:\공유 드라이브\BSL-Data\카이스트_단락셀\카이스트 단락셀\2차 셀 데이터\Aging\사이클 데이터';
files = dir(folder);

files = files(3:end,:);

%save in struct
data_st = struct('name',[],'cycle',[],'index',[]);
for i = 1:length(files)
    file_name = files(i).name;
    data_now = readtable([folder filesep file_name], "FileType","text","ReadVariableNames",false,'NumHeaderLines',11);

    %data trim (trim last cycle)
    data_now = data_now(data_now.Var15 ~= 0,:);

    %reassign name and data
    name_idx = strfind(file_name,'_');
    name_doe = file_name(1:name_idx(2) - 1);
    cycle_num = data_now.Var2;
    cycle_max = max(cycle_num);

    %save in struct
    data_st(i).name = name_doe;
    data_st(i).cycle = cycle_max;
end

%%Data separation upon DOE
%separate by DOE name
n = 1;

for k = 1:length(data_st)-1

    % index assign
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
        cycle_doe(f) = data_cal(f).cycle;
        name = split(data_cal(f).name,'_');
        name = name{1};
    end
    mean_data = mean(cycle_doe);
    std_data = std(cycle_doe);
    
    data_sum(j).name = name;
    data_sum(j).mean = mean_data;
    data_sum(j).std = std_data;
end

%plot
figure(1)
errorbar([data_sum.mean], [data_sum.std], 'bo', 'MarkerSize',6,'LineWidth',1.5, 'MarkerFaceColor','b')

%x-axis set
set(gca,'XTick',1:length({data_sum.name}));
set(gca,'XTickLabel',{data_sum.name});
xlim([0 length({data_sum.name})+1]);
ylim([0 22])

xtickangle(45)
xlabel('DOE');
ylabel('Max Cycle');

grid on;
box on;
