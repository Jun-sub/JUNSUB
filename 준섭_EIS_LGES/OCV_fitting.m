
clear, clc, close all; 


folder_path = 'G:\공유 드라이브\BSL-Data\LGES\2차 실험\OCP\데이터 변환'
filelist = dir([folder_path filesep '*.txt']);
filename = filelist(3)

for i = 1:length(filelist)

file_name = folder(i).name;
full = fullfile([folder_path filesep file_name])';

data = readtable(full, 'NumHeaderLines',10, 'ReadVariableNames',true, 'VariableNamingRule','preserve');


time = data(:,2);
current = data(:,7);
voltage = data(:,8);

cycle = data(:,3);

end 

