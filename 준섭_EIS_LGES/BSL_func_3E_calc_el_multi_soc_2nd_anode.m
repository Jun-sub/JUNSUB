function [f_data, z_integ_data, z_model, z_model_dist, paras_integ ,paras_integ_dist] = BSL_func_3E_calc_el_multi_soc_2nd_anode(path_folder, path_file, multi_soc_range, cell_type, anode_type)
% This code concurrently fits multiple SOC values to extract Del and Kel,
% minimizing Kel & Del discrepancies across different SOCs.

% NOTE: thie code accept [freq real imag] data not -imag 
% give anode_type only for 2nd anode fitting (not yet implemented)
% Func should locate before P2D fitting
% change name_file when you apply other files especially ay type_Acf 3 & 4
if length(multi_soc_range) == 1
    error('Should set more than 2 soc ranges')
end 
% data_load

type_acf = evalin("base",'type_acf'); %Call type_acf by evalin for preventing error with initial code


%% Anode & Cathode
if type_acf == 0 || type_acf == 1 || type_acf == 2
    name_file = cell(length(multi_soc_range),1);
    for i = 1:length(multi_soc_range)

         name_file{i} = sprintf(path_file, anode_type, multi_soc_range(i));
    
         data_raw = readmatrix([path_folder filesep name_file{i}]);
         data(:,3*i-2) = data_raw(:,1);
         if type_acf == 0
            data(:,3*i-1:3*i) = [data_raw(:,2) -data_raw(:,3)];
         elseif type_acf == 1
            data(:,3*i-1:3*i) = [-data_raw(:,4) data_raw(:,5)];
         elseif type_acf == 2
            data(:,3*i-1:3*i) = [data_raw(:,6) -data_raw(:,7)];
         end
    end
    
    if type_acf == 0
       data(18,:) = []; % trim full data in row 18. (failed data)
    end 
% data = [soc(1)freq soc(1)real soc(1)imag soc(2)freq soc(2)real soc(2)imag ... soc(i)freq soc(i)real soc(i)imag] 
% data identify and process upon each soc
 
    for i = 1:length(multi_soc_range)

        %trim inductance data and reconstruct 
        data_trim(:, 3*i-2:3*i) = data(data(:, 3*i) <= 0, 3*i-2:3*i);


        if i < length(multi_soc_range)
            if any(data(:,3*i-2) ~= data(:,3*(i+1)-2))
                error('%s soc and %s frequency data do not matched\nshould set narrower soc range',num2str(multi_soc_range(i)),num2str(multi_soc_range(i+1)))
            end 
        end

        %data_trim rearrangement (integrate freq data)
        z_integ_data(:,2*i-1:2*i) = data_trim(:,3*i-1:3*i);
    end

    f_data = data_trim(:,1);

   
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
    weight = zeros(length(z_integ_data(:,1)),2*length(multi_soc_range));

   if type_weight == 1  % relative weight
       
        for i = 1:length(multi_soc_range)
            weight(:, 2*i-1) = (z_integ_data(:, 2*i-1).^2 + z_integ_data(:, 2*i).^2).^(-0.5);
            weight(:, 2*i) = weight(:, 2*i-1);
        end
    
   elseif type_weight == 0  % absolute weight
        weight = ones(length(z_integ_data(:,1)),2*length(multi_soc_range));
   end

   fprintf('Start SOC integrated fitting \n')

%% Call EIS model
factors_integ_ini = ones(length(factors_ini)-2,length(multi_soc_range));
factors_integ_ini(1:2,end+1) = [1;1]; % for Kel & Del

lb = factors_integ_ini*0.0001;
ub = factors_integ_ini*1000;

weighted_model = @(factors,f_data)BSL_func_soc_integrated_model_V1_half(f_data,factors,multi_soc_range,T,type_acf,cell_type).*weight;
weighted_data = z_integ_data.*weight;

tic;
factors_integ_hat = lsqcurvefit(weighted_model,factors_integ_ini,f_data,weighted_data, lb, ub, options);
toc;

[z_model, paras_integ] = BSL_func_soc_integrated_model_V1_half(f_data,factors_integ_hat,multi_soc_range,T,type_acf,cell_type);

%% Call EIS + Dist model 
type_dist = evalin('base','type_dist');

factors_integ_ini_dist(1:length(factors_ini)-2,1:length(multi_soc_range)+1) = factors_integ_hat;
factors_integ_ini_dist(length(factors_ini)-1:length(factors_ini),1:length(multi_soc_range)) = 0.5; % set initial drt & ddt Std

weighted_model_dist = @(factors,f_data)BSL_func_EISmodel_V_half_soc_and_Dist_integrated(f_data,factors,multi_soc_range,T,type_acf,cell_type,type_dist).*weight;

lb = factors_integ_ini_dist*0.01;
ub = factors_integ_ini_dist*100;

fprintf('Start SOC integrated + Dist fitting \n')
tic;
   factors_integ_hat_dist = lsqcurvefit(weighted_model_dist,factors_integ_ini_dist,f_data,weighted_data, lb, ub, options_dist);
toc;

[z_model_dist, paras_integ_dist] = BSL_func_EISmodel_V_half_soc_and_Dist_integrated(f_data,factors_integ_hat_dist,multi_soc_range,T,type_acf,cell_type,type_dist);
% Check optimized data
    figure(length(multi_soc_range)*3) % for preventing of figure number overlap with main figure

t = tiledlayout(ceil(sqrt(length(multi_soc_range))),ceil(sqrt(length(multi_soc_range))),"TileSpacing","loose","Padding","loose");
for i = 1:length(multi_soc_range)
    nexttile 
    plot(z_integ_data(:,2*i-1),-z_integ_data(:,2*i), z_model(:,2*i-1),-z_model(:,2*i),z_model_dist(:,2*i-1),-z_model_dist(:,2*i))
    legend(['z data' ' soc ' num2str(multi_soc_range(i))],['z model' ' soc ' num2str(multi_soc_range(i))],['z model dist' ' soc ' num2str(multi_soc_range(i))])
    title(cell_type)
    axis equal
     set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')
end

set(gcf,'Position',[100 100 1200 1200]);

% Kel_factor_hat = factors_integ_hat_dist(1,end);
% Del_factor_hat = factors_integ_hat_dist(2,end);
% paras_integ_hat = paras_integ_dist; %Kel, Del

%% 3E_Sum
elseif type_acf == 3
    name_file = cell(2*length(multi_soc_range),1);
        for i = 1:2*length(multi_soc_range)
           
                if i <= length(multi_soc_range)
                    name_file{i} = sprintf(path_file, anode_type , multi_soc_range(i)); %Anode 
                    data_raw = readmatrix([path_folder filesep name_file{i}]);
                    data(:,3*i-2) = data_raw(:,1);
                    data(:,3*i-1:3*i) = [-data_raw(:,4) data_raw(:,5)];
                elseif i > length(multi_soc_range)
                    name_file{i} = sprintf(path_file, anode_type, multi_soc_range(i-length(multi_soc_range))); %Cathode data load
                    data(:,3*i-2) = data_raw(:,1);
                    data(:,3*i-1:3*i) = [data_raw(:,6) -data_raw(:,7)];
                end
                 
        end
         
        for i = 1:2*length(multi_soc_range)
            % Check freq data match
            if i < 2*length(multi_soc_range)
                if any(data(:,3*i-2) ~= data(:,3*(i+1)-2))
                    error('%s soc and %s frequency data do not matched\nshould set narrower soc range',num2str(multi_soc_range(i)),num2str(multi_soc_range(i+1)))
                end 
            end
            
            % Sum anode and cathode data
            if i <= length(multi_soc_range)
                z_data_sum(:,3*i-2) = data(:,3*i-2);
                z_data_sum(:,3*i-1:3*i) = data(:,3*i-1:3*i) + data(:,3*(i+length(multi_soc_range))-1:3*(i+length(multi_soc_range)));
                
                
                %trim inductance data and reconstruct
                data_trim(:, 3*i-2:3*i) = z_data_sum(data(:, 3*i) <= 0, 3*i-2:3*i);
                
        
                %data_trim rearrangement (integrate freq data)
                z_integ_data(:,2*i-1:2*i) = data_trim(:,3*i-1:3*i);
                % anode [real(i) imag(1) ... real(i) imag(i)]
            end
    
           
        end
    
            f_data = data_trim(:,1);
    
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
            weight = zeros(length(z_integ_data(:,1)),2*length(multi_soc_range));
        
           if type_weight == 1  % relative weight
               
                for i = 1:length(multi_soc_range)
                    weight(:, 2*i-1) = (z_integ_data(:, 2*i-1).^2 + z_integ_data(:, 2*i).^2).^(-0.5);
                    weight(:, 2*i) = weight(:, 2*i-1);
                end
            
           elseif type_weight == 0  % absolute weight
                weight = ones(length(z_integ_data(:,1)),2*length(multi_soc_range));
           end
        
           fprintf('Start SOC integrated fitting \n')
        
        %% Call EIS model
        factors_integ_ini = ones((length(factors_ini)-2)/2,2*length(multi_soc_range)); %factors will be apllied in integrated form
        factors_integ_ini(1:2,end+1) = [1;1]; % for Kel & Del
    
        
        lb = factors_integ_ini*0.0001;
        ub = factors_integ_ini*1000;
        
        weighted_model = @(factors,f_data)BSL_func_soc_integrated_model_V1_3E(f_data,factors,multi_soc_range,T,type_acf,cell_type).*weight;
        weighted_data = z_integ_data.*weight;
        
        tic;
        factors_integ_hat = lsqcurvefit(weighted_model,factors_integ_ini,f_data,weighted_data, lb, ub, options);
        toc;
        
        [z_model, paras_integ] = BSL_func_soc_integrated_model_V1_3E(f_data,factors_integ_hat,multi_soc_range,T,type_acf,cell_type);
        
        %% Call EIS + Dist model 
        type_dist = evalin('base','type_dist');
        
        factors_integ_ini_dist(1:(length(factors_ini)-2)/2,1:2*length(multi_soc_range)+1) = factors_integ_hat;
        factors_integ_ini_dist((length(factors_ini)-2)/2+1:(length(factors_ini)-2)/2+2,1:2*length(multi_soc_range)) = 0.5; % set initial drt & ddt Std
        
        weighted_model_dist = @(factors,f_data)BSL_func_EISmodel_V_3E_soc_and_Dist_integrated(f_data,factors,multi_soc_range,T,type_acf,cell_type,type_dist).*weight;
        
        lb = factors_integ_ini_dist*0.01;
        ub = factors_integ_ini_dist*100;
        
        fprintf('Start SOC integrated + Dist fitting \n')
        tic;
           factors_integ_hat_dist = lsqcurvefit(weighted_model_dist,factors_integ_ini_dist,f_data,weighted_data, lb, ub, options_dist);
        toc;
        
        [z_model_dist, paras_integ_dist] = BSL_func_EISmodel_V_3E_soc_and_Dist_integrated(f_data,factors_integ_hat_dist,multi_soc_range,T,type_acf,cell_type,type_dist);
        % Check optimized data
            figure(length(multi_soc_range)*3) % for preventing of figure number overlap with main figure
        
        t = tiledlayout(ceil(sqrt(length(multi_soc_range))),ceil(sqrt(length(multi_soc_range))),"TileSpacing","loose","Padding","loose");
        for i = 1:length(multi_soc_range)
            nexttile 
            plot(z_integ_data(:,2*i-1),-z_integ_data(:,2*i), z_model(:,2*i-1),-z_model(:,2*i),z_model_dist(:,2*i-1),-z_model_dist(:,2*i))
            legend(['z data' ' soc ' num2str(multi_soc_range(i))],['z model' ' soc ' num2str(multi_soc_range(i))],['z model dist' ' soc ' num2str(multi_soc_range(i))])
            title(cell_type)
            axis equal
             set(gca,'Box','on',... %Axis Properties: BOX   
            'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
            'FontUnits','points','FontSize',10,'FontName','Times New Roman')
            grid on;
            xlabel('Z_{re} [Ohm]')
            ylabel('-Z_{im} [Ohm]')
        end
        
        set(gcf,'Position',[100 100 1200 1200]);
        
        % Kel_factor_hat = factors_integ_hat_dist(1,end);
        % Del_factor_hat = factors_integ_hat_dist(2,end);
        % paras_integ_hat = paras_integ_dist; %Kel, Del
%% 3E_Simul
elseif type_acf == 4 
    name_file = cell(2*length(multi_soc_range),1);
        for i = 1:2*length(multi_soc_range)
           
                if i <= length(multi_soc_range)
                    name_file{i} = sprintf(path_file, anode_type , multi_soc_range(i)); %Anode 
                    data_raw = readmatrix([path_folder filesep name_file{i}]);
                    data(:,3*i-2) = data_raw(:,1);
                    data(:,3*i-1:3*i) = [-data_raw(:,4) data_raw(:,5)];
                elseif i > length(multi_soc_range)
                    name_file{i} = sprintf(path_file, anode_type, multi_soc_range(i-length(multi_soc_range))); %Cathode data load
                    data(:,3*i-2) = data_raw(:,1);
                    data(:,3*i-1:3*i) = [data_raw(:,6) -data_raw(:,7)];
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
            dataa_trim(:, 3*i-2:3*i) = data(data(:, 3*i) <= 0, 3*i-2:3*i);
        elseif i > length(multi_soc_range)
            datac_trim(:, 3*i-2:3*i) = data(data(:, 3*i) <= 0, 3*i-2:3*i);
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
    factors_integ_ini = ones((length(factors_ini)-2)/2,2*length(multi_soc_range)); %factors will be apllied in integrated form
    factors_integ_ini(1:2,end+1) = [1;1]; % for Kel & Del

    
    lb = factors_integ_ini*0.0001;
    ub = factors_integ_ini*1000;
    
    weighted_model = @(factors,f_data)BSL_func_soc_integrated_model_V1_3E(f_data,factors,multi_soc_range,T,type_acf,cell_type).*weight;
    weighted_data = z_integ_data.*weight;
    
    tic;
    factors_integ_hat = lsqcurvefit(weighted_model,factors_integ_ini,f_data,weighted_data, lb, ub, options);
    toc;
    
    [z_model, paras_integ] = BSL_func_soc_integrated_model_V1_3E(f_data,factors_integ_hat,multi_soc_range,T,type_acf,cell_type);
    % z_model = z_model_sep(:,1:end/2) + z_model_sep(:,end/2+1:end);
    %% Call EIS + Dist model 
    type_dist = evalin('base','type_dist');
    
    factors_integ_ini_dist(1:(length(factors_ini)-2)/2,1:2*length(multi_soc_range)+1) = factors_integ_hat;
    factors_integ_ini_dist((length(factors_ini)-2)/2+1:(length(factors_ini)-2)/2+2,1:2*length(multi_soc_range)) = 0.5; % set initial drt & ddt Std
    
    weighted_model_dist = @(factors,f_data)BSL_func_EISmodel_V_3E_soc_and_Dist_integrated(f_data,factors,multi_soc_range,T,type_acf,cell_type,type_dist).*weight;
    
    lb = factors_integ_ini_dist*0.01;
    ub = factors_integ_ini_dist*100;
    
    fprintf('Start SOC integrated + Dist fitting \n')
    tic;
       factors_integ_hat_dist = lsqcurvefit(weighted_model_dist,factors_integ_ini_dist,f_data,weighted_data, lb, ub, options_dist);
    toc;
    
    [z_model_dist, paras_integ_dist] = BSL_func_EISmodel_V_3E_soc_and_Dist_integrated(f_data,factors_integ_hat_dist,multi_soc_range,T,type_acf,cell_type,type_dist);
    % z_model_dist = z_model_dist_sep(:,1:end/2) + z_model_dist_sep(:,end/2+1:end);
    figure(length(multi_soc_range)*3) % for preventing of figure number overlap with main figure / Anode
        t = tiledlayout(ceil(sqrt(length(multi_soc_range))),ceil(sqrt(length(multi_soc_range))),"TileSpacing","loose","Padding","loose");
            for i = 1:length(multi_soc_range)
                nexttile 
                plot(z_integ_data(:,2*(i+length(multi_soc_range))-1),-z_integ_data(:,2*(i+length(multi_soc_range))), z_model(:,2*(i+length(multi_soc_range))-1),-z_model(:,2*(i+length(multi_soc_range))),z_model_dist(:,2*(i+length(multi_soc_range))-1),-z_model_dist(:,2*(i+length(multi_soc_range))))
                legend(['z data' ' soc ' num2str(multi_soc_range(i))],['z model' ' soc ' num2str(multi_soc_range(i))],['z model dist' ' soc ' num2str(multi_soc_range(i))])
                title('Anode')
                axis equal
                set(gca,'Box','on',... %Axis Properties: BOX   
                'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
                'FontUnits','points','FontSize',10,'FontName','Times New Roman')
                grid on;
                xlabel('Z_{re} [Ohm]')
                ylabel('-Z_{im} [Ohm]')
            end

    set(gcf,'Position',[100 100 1200 1200]);
    figure(length(multi_soc_range)*4) % for preventing of figure number overlap with main figure / Cathode
     t = tiledlayout(ceil(sqrt(length(multi_soc_range))),ceil(sqrt(length(multi_soc_range))),"TileSpacing","loose","Padding","loose");
            for i = 1:length(multi_soc_range)
                nexttile 
                plot(z_integ_data(:,2*(i+length(multi_soc_range))-1),-z_integ_data(:,2*(i+length(multi_soc_range))), z_model(:,2*(i+length(multi_soc_range))-1),-z_model(:,2*(i+length(multi_soc_range))),z_model_dist(:,2*(i+length(multi_soc_range))-1),-z_model_dist(:,2*(i+length(multi_soc_range))))
                legend(['z data' ' soc ' num2str(multi_soc_range(i))],['z model' ' soc ' num2str(multi_soc_range(i))],['z model dist' ' soc ' num2str(multi_soc_range(i))])
                title('Cathode')
                axis equal
                set(gca,'Box','on',... %Axis Properties: BOX   
                'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
                'FontUnits','points','FontSize',10,'FontName','Times New Roman')
                grid on;
                xlabel('Z_{re} [Ohm]')
                ylabel('-Z_{im} [Ohm]')
            end

    set(gcf,'Position',[1300 100 1200 1200]);

    % Kel_factor_hat = factors_integ_hat_dist(1,end);
    % Del_factor_hat = factors_integ_hat_dist(2,end);
    % paras_integ_hat = paras_integ_dist; %Kel, Del
end 
end 