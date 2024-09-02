%Aging visuallization

clear; clc; close all;

%데이터를 불러와서 이어서 플랏

% 폴더 지정하면 각 파일 분류할 수 있게.
full_folder_path = 'G:\공유 드라이브\BSL-Data\카이스트_단락셀\카이스트 단락셀\2차 셀 데이터\Aging'
full_folder = dir(full_folder_path);

num_files = length(full_folder);

overall_tic = tic;

for i = 3:1:num_files
fprintf('Processing file %d of %d: %s\n', i-2, num_files-2, full_folder(i).name);

step_tic = tic;

data_name = full_folder(i).name;
data_path = fullfile(full_folder_path, data_name);

data = readtable(data_path, "FileType","text",'NumHeaderLines', 11, 'ReadVariableNames', false);

% Data processing
t = data.Var2;
cycle = data.Var3;
I = data.Var7;
V = data.Var8;

figure(i-2)
yyaxis left
plot(t,I,'b')
ylabel ('Current [A]')
yyaxis right
plot(t,V,'r')
xlabel('Time [S]')
ylabel ('Voltage [V]')

name = strjoin(strsplit(data_name,'_'),' ');
title(name)

    step_toc = toc(step_tic);
    fprintf('Elapsed time for file %d: %.2f seconds\n', i-2, step_toc);
    
    overall_toc = toc(overall_tic);
    fprintf('Cumulative elapsed time: %.2f seconds\n\n', overall_toc);
% Cycle 별 데이터 필요할 때 사용
    % for j = 1:max(cycle)
    % 
    %  idx = (cycle == j);
    %  figure()
    %  yyaxis left
    %  plot(t(idx),I(idx),'b')
    %  ylabel ('Current [A]')
    %  yyaxis right
    %  plot(t(idx),V(idx),'r')
    %  xlabel('Time [S]')
    %  ylabel ('Voltage [V]')
    %  legend(['cycle', num2str(j)])
    % end
end

total_toc = toc(overall_tic);
    fprintf('Total elapsed time: %.2f seconds\n', total_toc);