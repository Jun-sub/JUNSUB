function [f_data, z_integ_data, z_model, z_model_dist, paras_integ ,paras_integ_dist] = BSL_func_3E_calc_el_multi_soc_V_final(path_folder, path_file, multi_soc_range, cell_type, anode_type)
% This code concurrently fits multiple SOC values to extract Del and Kel,
% minimizing Kel & Del discrepancies across different SOCs.

% NOTE: thie code accept [freq real imag] data not -imag 
% Func should locate before P2D fitting
% change name_file when you apply other files especially ay type_Acf 3 & 4
% if length(multi_soc_range) == 1
%     error('Should set more than 2 soc ranges')
% end 
% data_load

type_acf = evalin("base",'type_acf'); %Call type_acf by evalin for preventing error with initial code
%% 3E_Simul
    name_file = cell(2*length(multi_soc_range),1);
    for i = 1:2*length(multi_soc_range)
        if i <= length(multi_soc_range)
            name_file{i} = sprintf(path_file, 'Anode', multi_soc_range(i)); %Anode data load
            data(:,3*i-2:3*i) = readmatrix([path_folder filesep name_file{i}]);
        elseif i > length(multi_soc_range)
            name_file{i} = sprintf(path_file, 'Cathode', multi_soc_range(i-length(multi_soc_range))); %Catdhoe data loda
            data(:,3*i-2:3*i) = readmatrix([path_folder filesep name_file{i}]);
        end
    end
     
    for i = 1:2*length(multi_soc_range)

        if i < 2*length(multi_soc_range)
            if any(data(:,3*i-2) ~= data(:,3*(i+1)-2))
                error('%s soc and %s frequency data do not matched\nshould set narrower soc range',num2str(multi_soc_range(i)),num2str(multi_soc_range(i+1)))
            end 
        end
        
        %trim inductance data and reconstruct
        if i <= length(multi_soc_range)
            if i == 1
            dataa_trim(:, 3*i-2:3*i) = data(data(:, 3*i) <= 0, 3*i-2:3*i);
            end 

            if length(dataa_trim) ~= length(data(data(:, 3*i) <= 0, 3*i-2:3*i))
               leng_trim = length(dataa_trim);
               leng_data = length(data(:, 3*i));
               if leng_trim < leng_data
                  data_temp = data(leng_data-leng_trim+1:end,3*i-2:3*i);
                  dataa_trim(:, 3*i-2:3*i) = data_temp;
               elseif leng_trim > leng_data
                  dataa_trim(:, 3*i-2:3*i) = dataa_trim(end-leng_data+1:end,3*i-2:3*i);
               end 
            else 
               dataa_trim(:, 3*i-2:3*i) = data(data(:, 3*i) <= 0, 3*i-2:3*i);
            end
        elseif i > length(multi_soc_range)
            if i == 1 + length(multi_soc_range)
               datac_trim(:, 3*i-2:3*i) = data(data(:, 3*i) <= 0, 3*i-2:3*i);
            end

            if length(datac_trim) ~= length(data(data(:, 3*i) <= 0, 3*i-2:3*i))
               leng_trim = length(datac_trim);
               leng_data = length(data(:, 3*i));
               if leng_trim < leng_data
                  data_temp = data(leng_data-leng_trim+1:end,3*i-2:3*i);
                  datac_trim(:, 3*i-2:3*i) = data_temp;
               elseif leng_trim > leng_data
                  datac_trim(:, 3*i-2:3*i) = datac_trim(end-leng_data+1:end,3*i-2:3*i);
               end 
            else 
               datac_trim(:, 3*i-2:3*i) = data(data(:, 3*i) <= 0, 3*i-2:3*i);
            end
        end 

        %data_trim rearrangement (integrate freq data)
        if i <= length(multi_soc_range)
               z_integ_dataa(:,2*i-1:2*i) = dataa_trim(:,3*i-1:3*i);
               % anode [real(i) imag(1) ... real(i) imag(i)]
        elseif i > length(multi_soc_range)
            z_integ_datac(:,2*(i-length(multi_soc_range))-1:2*(i-length(multi_soc_range))) = datac_trim(:,3*i-1:3*i);
            % cathode [real(i) imag(1) ... real(i) imag(i)]
        end

       
    end

        % If the anode & cathode data length is not same -> equal the length
        if length(z_integ_datac(:,1)) ~= length(z_integ_dataa(:,1))
           
           leng_c = length(z_integ_datac(:,1));
           leng_a = length(z_integ_dataa(:,1));

           min_leng = min(leng_a,leng_c);
           min_start = abs(leng_a - leng_c)+1;

           % for making sure to have same freq range
           if leng_a < leng_c
               z_integ_dataa = z_integ_dataa(1:min_leng,:);
               z_integ_datac = z_integ_datac(min_start:end,:);
               f_data = dataa_trim(:,1);
           elseif leng_a > leng_c
               z_integ_dataa = z_integ_dataa(min_start:end,:);
               z_integ_datac = z_integ_datac(1:min_leng,:);
               f_data = datac_trim(:,3*length(multi_soc_range) + 1);
           end 
        end

        z_integ_data = [z_integ_dataa z_integ_datac];
        % z_integ_data = z_integ_data_sep(:,1:end/2) + z_integ_data_sep(:,end/2+1:end);
        

    %% Optimization 
    
    % call subtle elements from base work spaceb by using evalin
    % check all the elements are exist
    T = evalin("base",'T');
    
    type_weight = evalin("base",'type_weight');
    
    num_iter = evalin('base','num_iter');
    num_iter_dist = evalin('base','num_iter_dist');
    
    options = evalin("base",'options');
    options_dist = evalin("base",'options_dist');
    
    
    factors_ini = evalin("base",'factors_ini');
    
    
        if ~exist("T",'var') || ~exist('type_weight','var') || ~exist('options','var') || ~exist('options_dist', 'var') ||  ~exist('factors_ini', 'var') ||  ~exist('type_acf', 'var')
            error('Assign variables to function')
        end 
    
    % Define weighting vector
    % initial matrix set    
    weight = zeros(length(z_integ_data(:,1)),4*length(multi_soc_range));
    
       if type_weight == 1  % relative weight
           
            for i = 1:2*length(multi_soc_range)
                if i <= length(multi_soc_range)
                weight(:, 2*i-1) = (z_integ_data(:, 2*i-1).^2 + z_integ_data(:, 2*i).^2).^(-0.5);
                weight(:, 2*i) = weight(:, 2*i-1);
                elseif i > length(multi_soc_range)
                weight(:, 2*i-1) = (z_integ_data(:, 2*i-1).^2 + z_integ_data(:, 2*i).^2).^(-0.5);
                weight(:, 2*i) = weight(:, 2*i-1);
                end
            end
        
       elseif type_weight == 0  % absolute weight
            weight = ones(length(z_integ_dataa(:,1)),4*length(multi_soc_range));
       end
    
       fprintf('Start SOC integrated fitting \n')
    
    %% Call EIS model
    factors_integ_ini = ones((length(factors_ini)-2)/2,2*length(multi_soc_range)); %factors will be apllied in integrated form, 
    %1 = R_itsc; 2 = i0; 3 = Cdl; 4 = Ds; 5 = Av; 6 = DRT_std; 7 = DDT_std,
    %temp 8 = Kel; 9 = Del;
    factors_integ_ini(1:6,2*length(multi_soc_range)+1) = [1;1;1;1;4.098543764913605e-06;4.098543764913605e-06]; % for Kela; Dela; Av_n; Av_p; y_shift1(0-0.6); y_shift2(0.61-1)
    factors_integ_ini(1:6,2*length(multi_soc_range)+2) = [1;1;0.237;0.234;1;1]; %for Kelc; Delc;Epsla;Epslc;a_n5;a_p

    % temporary for test
    % load("C:\Users\admin\Documents\GitHub\JunSub\준섭_EIS_LGES\12월 미팅\dUdc_adjust_test\workspace_natural.mat")
    % factors_integ_ini = factors_integ_hat;
    % 
    % options= optimset('Display','iter','MaxIter',10,'MaxFunEvals',1e6,...
    %     'TolFun',1e-8,'TolX',1e-8,'FinDiffType','central');
    % ------------------

    lb = factors_integ_ini*0.001;
    ub = factors_integ_ini*200;
    
    % Linear gradient lb & ub 
    lb(1:6,2*length(multi_soc_range)+2) = 0; % for gradient
    ub(1:6,2*length(multi_soc_range)+2) = 30;
    
    % Kel & Del lb & ub
    lb(1:2,2*length(multi_soc_range)+1) = 0.01; % for Kel, Del
    ub(1:2,2*length(multi_soc_range)+1) = 20;
    
    weighted_model = @(factors,f_data)BSL_func_soc_integrated_model_V_final_3E(f_data,factors,multi_soc_range,T,type_acf,cell_type).*weight;
    weighted_data = z_integ_data.*weight;
    
    tic;
    factors_integ_hat = lsqcurvefit(weighted_model,factors_integ_ini,f_data,weighted_data, lb, ub, options);
    toc;
    
    [z_model, paras_integ] = BSL_func_soc_integrated_model_V_final_3E(f_data,factors_integ_hat,multi_soc_range,T,type_acf,cell_type);

    %% Call EIS + Dist model 
    type_dist = evalin('base','type_dist');
    
    factors_integ_ini_dist(1:6,1:2*length(multi_soc_range)+2) = factors_integ_hat;
    factors_integ_ini_dist(6:7,1:2*length(multi_soc_range)) = 0.5; % set initial drt & ddt Std
    
    weighted_model_dist = @(factors,f_data)BSL_func_EISmodel_V_final_3E_soc_and_Dist_integrated(f_data,factors,multi_soc_range,T,type_acf,cell_type,type_dist).*weight;
    
    %lb & ub set
    lb = factors_integ_ini_dist*0.1;
    ub = factors_integ_ini_dist*10;
    
    % % Eplsa & Epslc lb & ub
    % lb(5,1:2*length(multi_soc_range)) = 0.05; % for epsla
    % ub(5,1:2*length(multi_soc_range)) = 0.4;

    % DRT & DDT lb & ub 
    lb(6:7,1:2*length(multi_soc_range)) = 0.01;
    ub(6:7,1:2*length(multi_soc_range)) = 5;

    lb(5,1:2*length(multi_soc_range)) = 0; % for gradient
    ub(5,1:2*length(multi_soc_range)) = 30;


    fprintf('Start SOC integrated + Dist fitting \n')
    % tic;
    %    factors_integ_hat_dist = lsqcurvefit(weighted_model_dist,factors_integ_ini_dist,f_data,weighted_data, lb, ub, options_dist);
    % toc;
    
    %-----tentative
    factors_integ_hat_dist = factors_integ_ini_dist;

    [z_model_dist, paras_integ_dist] = BSL_func_EISmodel_V_final_3E_soc_and_Dist_integrated(f_data,factors_integ_hat_dist,multi_soc_range,T,type_acf,cell_type,type_dist);
    % z_model_dist = z_model_dist_sep(:,1:end/2) + z_model_dist_sep(:,end/2+1:end);
    

    assignin("base", "f_data","f_data")
    assignin("base","factors_integ_hat",factors_integ_hat)
    assignin("base","factors_integ_hat_dist",factors_integ_hat_dist)

    % Kel_factor_hat = factors_integ_hat_dist(1,end);
    % Del_factor_hat = factors_integ_hat_dist(2,end);
    % paras_integ_hat = paras_integ_dist; %Kel, Del
end 