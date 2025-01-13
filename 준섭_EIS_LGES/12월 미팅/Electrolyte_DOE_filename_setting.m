clc; clear; close all;

folder_path = 'G:\공유 드라이브\BSL-Data\LGES\전해질 DOE\전해액 변경_DOE 공유\processed_data';
folder = dir(folder_path);

folder = folder(3:end,:);

for i = 1:length(folder)

    new_name{i} = strrep(folder(i).name,'%','');
    

    file_old = fullfile(folder_path,folder(i).name);
    file_new = fullfile(folder_path,new_name{i});

    movefile(file_old,file_new)

end