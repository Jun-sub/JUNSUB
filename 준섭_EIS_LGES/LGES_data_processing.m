% This code is for processing raw EIS data from LGES anode DOE 2
% column 1: freq, 2-3: full, 4-5: anode, 6-7: cathode

% freq/Hz : 주파수
% 
% Re(Z)/Ohm : fullcell real
% 
% -Im(Z)/Ohm : fullcell imginary
% 
% Re(Zce)/Ohm: anode real
% 
% -Im(Zce)/Ohm: anode imginary
% 
% Re(Zwe-ce)/Ohm: cathode real

% target format (freq/ +real/ +imag)
clc, clear, close all;

% Data load
folder_path = 'G:\공유 드라이브\BSL-Data\LGES\전해질 DOE\전해액 변경_DOE 공유';
save_path = 'G:\공유 드라이브\BSL-Data\LGES\전해질 DOE\전해액 변경_DOE 공유\processed_data';

folder = dir(folder_path);
if folder(1).bytes == 0 
    folder = folder(4:end,1:end); %exclude processed_data folder
else 
    folder = folder(3:end,1:end); %exclude processed_data folder
end

for i = 1:length(folder)
    file_name = fullfile(folder_path, folder(i).name);
    
    data_raw = readmatrix(file_name);
    f_data = data_raw(:,1);

    for k = 1:3 %for separating full, anode, cathode
        
        data_save(:,1) = f_data;
        if data_raw(1,2*k) >= 0 %for making freq/real/imag format
           data_save(:,2) = data_raw(:,2*k);
        else
           data_save(:,2) = -data_raw(:,2*k);
        end 

        if data_raw(1,1+2*k) >= 0 %for making freq/real/imag format
           data_save(:,3) = data_raw(:,1+2*k); 
        else
           data_save(:,3) = -data_raw(:,1+2*k);
        end 
    
    
        save_name = split(folder(i).name,'.');
        
        if k == 1
           save_name = ['full_' save_name{1}];
        elseif k == 2
           save_name = ['anode_' save_name{1}];
        elseif k == 3
            save_name = ['cathode_' save_name{1}];
        end 

        save_path_full = fullfile(save_path,[save_name, '.csv']);
     
    
        writematrix(data_save,save_path_full)
    end 
end