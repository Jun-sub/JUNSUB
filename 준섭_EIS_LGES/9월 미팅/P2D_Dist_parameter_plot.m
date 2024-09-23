%This code stands for plot of parameters vs SOC 
%from BSL_Fitting_V1_Dist fitting

clc, clear, close all;

%data load
% set file folder upon your intention before running. 
% Ctrl + H --> DRT or DDT convert to each other.
load('C:\Users\admin\Documents\GitHub\JunSub\준섭_EIS_LGES\Dist_model\Cathode_DRT\para_sum.mat');



%data separation
soc = para_sum(1,:); %soc. common use
i0_P2D = para_sum(3,:); 
Ds_P2D = para_sum(5,:);
K_el_P2D = para_sum(6,:);
D_el_P2D = para_sum(7,:);

i0_dist = para_sum(12,:); 
Ds_dist = para_sum(14,:);
K_el_dist = para_sum(15,:);
D_el_dist = para_sum(16,:);

%%plot
%i0
figure(1)
t = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile
plot(soc, i0_P2D,'ro-','LineWidth',1,'MarkerFaceColor','r'); %Sum_anode
    set(gca,'Box','on',...    
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    title('P2D i0')
    xlabel('soc')
    ylabel('i0')
    legend('i0 P2D');

nexttile
plot(soc, i0_dist,'bo-','LineWidth',1,'MarkerFaceColor','b'); %Sum_cathode
set(gca,'Box','on',...  
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
     title('P2D + DRT i0')
    xlabel('soc')
    ylabel('i0')
    legend('i0 P2D + DRT');

nexttile
plot(soc, i0_P2D,'ro-','LineWidth',1,'MarkerFaceColor','r'); hold on
plot(soc, i0_dist,'bo-','LineWidth',1,'MarkerFaceColor','b'); %Sum_cathode
set(gca,'Box','on',...  
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    hold off
     title('P2D vs P2D + DRT')
    xlabel('soc')
    ylabel('i0')
    legend('i0 P2D','i0 P2D + DRT');

    set(gcf,'Position',[100 100 400 1200])

%Ds
figure(2)
t = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile
plot(soc, log10(Ds_P2D),'ro-','LineWidth',1,'MarkerFaceColor','r'); %Sum_anode
    set(gca,'Box','on',...    
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    title('P2D Ds')
    xlabel('soc')
    ylabel('log[Ds]')
    legend('Ds P2D');

nexttile
plot(soc, log10(Ds_dist),'bo-','LineWidth',1,'MarkerFaceColor','b'); %Sum_cathode
set(gca,'Box','on',...  
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
     title('P2D + DRT Ds')
    xlabel('soc')
    ylabel('log[Ds]')
    legend('Ds P2D + DRT');

nexttile
plot(soc, log10(Ds_P2D),'ro-','LineWidth',1,'MarkerFaceColor','r'); hold on
plot(soc, log10(Ds_dist),'bo-','LineWidth',1,'MarkerFaceColor','b'); %Sum_cathode
set(gca,'Box','on',...  
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    hold off
     title('P2D vs P2D + DRT')
    xlabel('soc')
    ylabel('log[Ds]')
    legend('Ds P2D','Ds P2D + DRT');

    set(gcf,'Position',[100 100 400 1200])

%K_el
figure(3)
t = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile
plot(soc, K_el_P2D,'ro-','LineWidth',1,'MarkerFaceColor','r'); %Sum_anode
    set(gca,'Box','on',...    
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    title('P2D Kel')
    xlabel('soc')
    ylabel('Kel')
    legend('Kel P2D');

nexttile
plot(soc, K_el_dist,'bo-','LineWidth',1,'MarkerFaceColor','b'); %Sum_cathode
set(gca,'Box','on',...  
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
     title('P2D + DRT Kel')
    xlabel('soc')
    ylabel('Kel')
    legend('Kel P2D + DRT');

nexttile
plot(soc, K_el_P2D,'ro-','LineWidth',1,'MarkerFaceColor','r'); hold on
plot(soc, K_el_dist,'bo-','LineWidth',1,'MarkerFaceColor','b'); %Sum_cathode
set(gca,'Box','on',...  
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    hold off
     title('P2D vs P2D + DRT')
    xlabel('soc')
    ylabel('Kel')
    legend('Kel P2D','Kel P2D + DRT');

    set(gcf,'Position',[100 100 400 1200])

    %D_el
figure(4)
t = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile
plot(soc, log10(D_el_P2D),'ro-','LineWidth',1,'MarkerFaceColor','r'); %Sum_anode
    set(gca,'Box','on',...    
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    title('P2D Del')
    xlabel('soc')
    ylabel('log[Del]')
    legend('Del P2D');

nexttile
plot(soc, log10(D_el_dist),'bo-','LineWidth',1,'MarkerFaceColor','b'); %Sum_cathode
set(gca,'Box','on',...  
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
     title('P2D + DRT Del')
    xlabel('soc')
    ylabel('log[Del]')
    legend('Del P2D + DRT');

nexttile
plot(soc, log10(D_el_P2D),'ro-','LineWidth',1,'MarkerFaceColor','r'); hold on
plot(soc, log10(D_el_dist),'bo-','LineWidth',1,'MarkerFaceColor','b'); %Sum_cathode
set(gca,'Box','on',...  
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    hold off
     title('P2D vs P2D + DRT')
    xlabel('soc')
    ylabel('log[Del]')
    legend('Del P2D','Del P2D + DRT');

    set(gcf,'Position',[100 100 400 1200])