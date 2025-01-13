%This code is for separating the size of each component (R_ct, (real)R_diff,
%(imag)R_diff

clc, clear, close all;

% Configuration
path_1 = 'G:\공유 드라이브\BSL-Data\LGES\LG raw data\12_6cm2_soc10_EIS # Sample 1'; % base
path_2 = 'G:\공유 드라이브\BSL-Data\LGES\Modeling data_2차 음극 2종\processed_data'; % Blend
path_3 = 'G:\공유 드라이브\BSL-Data\LGES\Modeling data_2차 음극 2종\processed_data'; % Natural

name_1 = 'PEIS_C09_%s_cycle_soc%d.csv'; % base
name_2 = '인조_천연 Blending음극_#1_soc%d_%s.csv'; % blend
name_3 = '천역흑연100_음극_#2_soc%d_%s.csv'; % natural

% save_path = 

type_acf = 1; % 1:anode, 2:cathode
soc = 20:10:80;
name_anode = {'Base' 'Blend' 'Natural'};
name_z = {'R_ct' 'R_diff' 'R_high_imag'};

if type_acf == 1
   cell_type = 'anode';
elseif type_acf == 2
   cell_type = 'cathode';
end

% load data
for i = 1:length(soc)
    data_base(:,3*i-2:3*i) = readmatrix(fullfile(path_1,sprintf(name_1,cell_type,soc(i)))); %base data
    data_blend(:,3*i-2:3*i) = readmatrix(fullfile(path_2,sprintf(name_2,soc(i),cell_type))); %blend data
    data_natural(:,3*i-2:3*i) = readmatrix(fullfile(path_3,sprintf(name_3,soc(i),cell_type))); %natural data
end

% Calculation of each elements
% 1) R_ct
% R_real(w_transition) - R_real(w_highest) 
% 2) real(R_diff)
% R_real(w_lowest) - R_real(w_transition)
% 3) imag(R_diff)
% R_imag(w_lowest)
for i = 1:length(soc)
    R_tran_base(i) = data_base(data_base(:,3*i) == max(data_base(data_base(:,3*i-2) < 1000,3*i)),3*i-1);
    R_tran_blend(i) = data_blend(data_blend(:,3*i) == max(data_blend(data_blend(:,3*i-2) < 1000,3*i)),3*i-1);
    R_tran_natural(i) = data_natural(data_natural(:,3*i) == max(data_natural(data_natural(:,3*i-2) < 1000,3*i)),3*i-1);

    R_high_base(i) = data_base(find(data_base(:,3*i) < 0,1),3*i-1);
    R_high_blend(i) = data_blend(find(data_blend(:,3*i) < 0,1),3*i-1);
    R_high_natural(i) = data_natural(find(data_natural(:,3*i) < 0,1),3*i-1);

    R_lower_r_base(i) = max(data_base(1:end-3,3*i-1));
    R_lower_r_blend(i) = max(data_blend(1:end-3,3*i-1));
    R_lower_r_natural(i) = max(data_natural(1:end-3,3*i-1));

    R_lower_i_base(i) = abs(data_base(end,3*i));
    R_lower_i_blend(i) = abs(data_blend(end,3*i));
    R_lower_i_natural(i) = abs(data_natural(end,3*i));
    
    % Calculation
    R_ct_base(i) = R_tran_base(i) - R_high_base(i);
    R_ct_blend(i) = R_tran_blend(i) - R_high_blend(i);
    R_ct_natural(i) = R_tran_natural(i) - R_high_natural(i);

    R_diff_base(i) = R_lower_r_base(i) - R_tran_base(i);
    R_diff_blend(i) = R_lower_r_blend(i) - R_tran_blend(i);
    R_diff_natural(i) = R_lower_r_natural(i) - R_tran_natural(i);
end

% dUdC calculation
for k = 1:length(soc)
    x_1 = 0.8781; % anode stoic when soc =0
    x_0 = 0.0216; % anode stoic when soc =1
    y_0 = 0.9319; % cathode stoic when soc =0
    y_1 = 0.3532; % cathode stoic when soc =1
    
    x = x_0 + (x_1 - x_0)*0.01*soc(k); % anode stoic
    y = y_0 + (y_1 - y_0)*0.01*soc(k); % cathode stoic

    cta = 29626;                        ctc = 48786;            % [mol/m3]

    dx = 0.0001; % finite difference step size.
    chg = 0.5; % amount weighting on charging curve wrpt discharging.
    dUdcc(k) = (1/ctc)*(Uc_function_v2(y+dx,chg) - Uc_function_v2(y-dx,chg))/(2*dx);   % {modified}
    dUdca(k) = (1/cta)*(Ua_function_v2(x+dx,chg) - Ua_function_v2(x-dx,chg))/(2*dx);    % *+*
      
    plateau1_start = 0.07; %for dUdca plateau setting
    plateau1_end = 0.15;

    plateau2_start = 0.18;
    plateau2_end = 0.50;

    plateau3_start = 0.51;
    plateau3_end = 0.98;



    if (plateau1_start <= soc(k)*0.01) && (soc(k)*0.01 <= plateau1_end) %soc 0.07 - 0.15
       dUdca(k) = (1/cta)*(Ua_function_v2(plateau1_end,chg) - Ua_function_v2(plateau1_start,chg))/(plateau1_end-plateau1_start);
    elseif (plateau2_start <= soc(k)*0.01) && (soc(k)*0.01 <= plateau2_end) % soc 0.18 - 0.50
       dUdca(k) = (1/cta)*(Ua_function_v2(plateau2_end,chg) - Ua_function_v2(plateau2_start,chg))/(plateau2_end-plateau2_start);
    elseif (plateau3_start <= soc(k)*0.01) && (soc(k)*0.01 <= plateau3_end) % 0.51 - 00.98
       dUdca(k) = (1/cta)*(Ua_function_v2(plateau3_end,chg) - Ua_function_v2(plateau3_start,chg))/(plateau3_end-plateau3_start);
    end 

    if type_acf == 1
        dUdc(k) = dUdca(k);
    elseif type_acf == 2
        dUdc(k) = dUdcc(k);
    end
end


%plot the value of each component with bar
figure()
t = tiledlayout(3,3,"TileSpacing","compact","Padding","compact");
%R_ct
    nexttile
    bar(soc,R_ct_base)
    set(gca,'PlotBoxAspectRatio',[1 1 1],'FontUnits','points','FontSize',10,'FontName','TimesNewRoman','Box','on')
    xlabel('soc')
    ylabel("|Z'|")
    ylim([0 0.35])
    legend(name_z{1})
    title('Base')
    
    nexttile
    bar(soc,R_ct_blend)
    set(gca,'PlotBoxAspectRatio',[1 1 1],'FontUnits','points','FontSize',10,'FontName','TimesNewRoman','Box','on')
    xlabel('soc')
    ylabel("|Z'|")
    ylim([0 0.35])
    legend(name_z{1})
    title('Blend')

    nexttile
    bar(soc,R_ct_natural)
    set(gca,'PlotBoxAspectRatio',[1 1 1],'FontUnits','points','FontSize',10,'FontName','TimesNewRoman','Box','on')
    xlabel('soc')
    ylabel("|Z'|")
    ylim([0 0.35])
    legend(name_z{1})
    title('natural')
% R_diff
    nexttile
    bar(soc,R_diff_base,"FaceColor",'r')
    set(gca,'PlotBoxAspectRatio',[1 1 1],'FontUnits','points','FontSize',10,'FontName','TimesNewRoman','Box','on')
    xlabel('soc')
    ylabel("|Z'|")
    ylim([0 1.85])
    legend(name_z{2})
    title('Base')

    nexttile
    bar(soc,R_diff_blend,"FaceColor",'r')
    set(gca,'PlotBoxAspectRatio',[1 1 1],'FontUnits','points','FontSize',10,'FontName','TimesNewRoman','Box','on')
    xlabel('soc')
    ylabel("|Z'|")
    ylim([0 1.85])
    legend(name_z{2})
    title('Blend')

    nexttile
    bar(soc,R_diff_natural,"FaceColor",'r')
    set(gca,'PlotBoxAspectRatio',[1 1 1],'FontUnits','points','FontSize',10,'FontName','TimesNewRoman','Box','on')
    xlabel('soc')
    ylabel("|Z'|")
    ylim([0 1.85])
    legend(name_z{2})
    title('natural')

% R_high_imag
    nexttile
    bar(soc,R_lower_i_base,"FaceColor",'y')
    set(gca,'PlotBoxAspectRatio',[1 1 1],'FontUnits','points','FontSize',10,'FontName','TimesNewRoman','Box','on')
    xlabel('soc')
    ylabel('|Z"|')
    ylim([0 8])
    legend(name_z{3})
    title('Base')

    nexttile
    bar(soc,R_lower_i_blend,"FaceColor",'y')
    set(gca,'PlotBoxAspectRatio',[1 1 1],'FontUnits','points','FontSize',10,'FontName','TimesNewRoman','Box','on')
    xlabel('soc')
    ylabel('|Z"|')
    ylim([0 8])
    legend(name_z{3})
    title('Blend')

    nexttile
    bar(soc,R_lower_i_natural,"FaceColor",'y')
    set(gca,'PlotBoxAspectRatio',[1 1 1],'FontUnits','points','FontSize',10,'FontName','TimesNewRoman','Box','on')
    xlabel('soc')
    ylabel('|Z"|')
    ylim([0 8])
    legend(name_z{3})
    title('natural')

    set(gcf,"Position",[300 100 1200 1200])

% R_D/dUdc plot
figure()
t = tiledlayout(2,3,"TileSpacing","compact","Padding","compact");

% R_diff
    nexttile
    bar(soc,-R_diff_base./dUdc,"FaceColor",'r')
    set(gca,'PlotBoxAspectRatio',[1 1 1],'FontUnits','points','FontSize',10,'FontName','TimesNewRoman','Box','on')
    xlabel('soc')
    ylabel("|Z'|")
    ylim([0 3.5e5])
    legend([name_z{2} '/dUdc'])
    title('Base')

    nexttile
    bar(soc,-R_diff_blend./dUdc,"FaceColor",'r')
    set(gca,'PlotBoxAspectRatio',[1 1 1],'FontUnits','points','FontSize',10,'FontName','TimesNewRoman','Box','on')
    xlabel('soc')
    ylabel("|Z'|")
    ylim([0 3.5e5])
    legend([name_z{2} '/dUdc'])
    title('Blend')

    nexttile
    bar(soc,-R_diff_natural./dUdc,"FaceColor",'r')
    set(gca,'PlotBoxAspectRatio',[1 1 1],'FontUnits','points','FontSize',10,'FontName','TimesNewRoman','Box','on')
    xlabel('soc')
    ylabel("|Z'|")
    ylim([0 3.5e5])
    legend([name_z{2} '/dUdc'])
    title('natural')

% R_high_imag
    nexttile
    bar(soc,-R_lower_i_base./dUdc,"FaceColor",'y')
    set(gca,'PlotBoxAspectRatio',[1 1 1],'FontUnits','points','FontSize',10,'FontName','TimesNewRoman','Box','on')
    xlabel('soc')
    ylabel('|Z"|')
    ylim([0 2e6])
    legend([name_z{3} '/dUdc'])
    title('Base')

    nexttile
    bar(soc,-R_lower_i_blend./dUdc,"FaceColor",'y')
    set(gca,'PlotBoxAspectRatio',[1 1 1],'FontUnits','points','FontSize',10,'FontName','TimesNewRoman','Box','on')
    xlabel('soc')
    ylabel('|Z"|')
    ylim([0 2e6])
    legend([name_z{3} '/dUdc'])
    title('Blend')

    nexttile
    bar(soc,-R_lower_i_natural./dUdc,"FaceColor",'y')
    set(gca,'PlotBoxAspectRatio',[1 1 1],'FontUnits','points','FontSize',10,'FontName','TimesNewRoman','Box','on')
    xlabel('soc')
    ylabel('|Z"|')
    ylim([0 2e6])
    legend([name_z{3} '/dUdc'])
    title('natural')

    set(gcf,"Position",[300 100 1200 800])

