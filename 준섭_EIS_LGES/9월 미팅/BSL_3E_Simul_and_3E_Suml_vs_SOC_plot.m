%This code stands for plot of parameters vs SOC 
%from 3E_simul & 3E_Sum fitting

%% 기존 3E_sum --> 3E_simul / 기존 3E2full --> 3E_sum으로 변경 필요
clc, clear, close all;

%data load

load('C:\Users\admin\Documents\GitHub\JunSub\준섭_EIS_LGES\3E2Sum & 3ESimul\para_sum3E_Simul.mat');
load('C:\Users\admin\Documents\GitHub\JunSub\준섭_EIS_LGES\3E2Sum & 3ESimul\para_sum3E_sum.mat');

data_3E_sum = para_sum3E_sum;
data_3E2full = para_sum3E2full;

%%data separation
%3E_sum
soc = data_3E_sum(1,:); %soc. common use

i0_sum = [data_3E_sum(3,:);data_3E_sum(8,:)]; %anode & cathode
Ds_sum = [data_3E_sum(5,:);data_3E_sum(10,:)]; %anode & cathode
av_sum = [data_3E_sum(6,:);data_3E_sum(11,:)];
K_el_sum = [data_3E_sum(12,:)]; %electrolyte
D_el_sum = [data_3E_sum(13,:)]; %electrolyte

%3E2full_sum
i0_3E2full = [data_3E2full(3,:);data_3E2full(8,:)]; %anode & cathode
Ds_3E2full = [data_3E2full(5,:);data_3E2full(10,:)]; %anode & cathode
av_3E2full = [data_3E2full(6,:);data_3E2full(11,:)];
K_el_3E2full = [data_3E2full(12,:)]; %electrolyte
D_el_3E2full = [data_3E2full(13,:)]; %electrolyte


%%plot
%av
figure(5)
t = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

nexttile
plot(soc, av_sum(1,:),'bo-','LineWidth',1,'MarkerFaceColor','b'); %Sum_anode
    set(gca,'Box','on',...    
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    title('3E_Simul anode')
    xlabel('soc')
    ylabel('av')
    legend('av Simul anode');

nexttile
plot(soc, av_sum(2,:),'bo-','LineWidth',1,'MarkerFaceColor','b'); %Sum_cathode
set(gca,'Box','on',...  
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
     title('3E simul cathode')
    xlabel('soc')
    ylabel('av')
    legend('av simul cathode');

nexttile
plot(soc, av_3E2full(1,:),'ro-','LineWidth',1,'MarkerFaceColor','r'); %3E2full_anode
set(gca,'Box','on',... 
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    hold off
    title('3E_Sum anode')
    xlabel('soc')
    ylabel('av')
    legend('av Sum anode');

nexttile
plot(soc, av_3E2full(2,:),'ro-','LineWidth',1,'MarkerFaceColor','r'); %3E2full_cathode
set(gca,'Box','on',...   
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
   title('3E Sum cathode')
    xlabel('soc')
    ylabel('av')
    legend('av 3E Sum cathode');

nexttile
plot(soc, av_sum(1,:),'bo-','LineWidth',1,'MarkerFaceColor','b'); hold on %Sum + 3E2full anode in one fig
plot(soc, av_3E2full(1,:),'ro-','LineWidth',1,'MarkerFaceColor','r'); %3E2full_anode
    set(gca,'Box','on',...    
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    hold off
    title('3E Simul + 3E Sum av anode')
    xlabel('soc')
    ylabel('av')
    legend('av Simul anode', 'av 3E Sum anode');

    nexttile
plot(soc, av_sum(2,:),'bo-','LineWidth',1,'MarkerFaceColor','b'); hold on%Sum + 3E2full cathode in one fig
plot(soc, av_3E2full(2,:),'ro-','LineWidth',1,'MarkerFaceColor','r');
set(gca,'Box','on',...  
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    hold off
        title('3E Simul + 3E Sum av cathode')
    xlabel('soc')
    ylabel('av')
    legend('av Simul cathode', 'av 3E Sum cathode');

set(gcf,"Position",[100 100 800 1200])


%------------




%%plot
%i0
figure(1)
t = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

nexttile
plot(soc, i0_sum(1,:),'bo-','LineWidth',1,'MarkerFaceColor','b'); %Sum_anode
    set(gca,'Box','on',...    
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    title('3Esum anode')
    xlabel('soc')
    ylabel('i0')
    legend('i0 sum anode');

nexttile
plot(soc, i0_sum(2,:),'bo-','LineWidth',1,'MarkerFaceColor','b'); %Sum_cathode
set(gca,'Box','on',...  
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
     title('3Esum cathode')
    xlabel('soc')
    ylabel('i0')
    legend('i0 sum cathode');

nexttile
plot(soc, i0_3E2full(1,:),'ro-','LineWidth',1,'MarkerFaceColor','r'); %3E2full_anode
set(gca,'Box','on',... 
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    hold off
    title('3E2full anode')
    xlabel('soc')
    ylabel('i0')
    legend('i0 3E2full anode');

nexttile
plot(soc, i0_3E2full(2,:),'ro-','LineWidth',1,'MarkerFaceColor','r'); %3E2full_cathode
set(gca,'Box','on',...   
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
   title('3E2full cathode')
    xlabel('soc')
    ylabel('i0')
    legend('i0 3E2full cathode');

nexttile
plot(soc, i0_sum(1,:),'bo-','LineWidth',1,'MarkerFaceColor','b'); hold on %Sum + 3E2full anode in one fig
plot(soc, i0_3E2full(1,:),'ro-','LineWidth',1,'MarkerFaceColor','r'); %3E2full_anode
    set(gca,'Box','on',...    
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    hold off
    title('3Esum + 3E2full i0 anode')
    xlabel('soc')
    ylabel('i0')
    legend('i0 sum anode', 'i0 3E2full anode');

    nexttile
plot(soc, i0_sum(2,:),'bo-','LineWidth',1,'MarkerFaceColor','b'); hold on%Sum + 3E2full cathode in one fig
plot(soc, i0_3E2full(2,:),'ro-','LineWidth',1,'MarkerFaceColor','r');
set(gca,'Box','on',...  
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    hold off
     title('3Esum + 3E2full i0 cathode')
    xlabel('soc')
    ylabel('i0') 
    legend('i0 sum cathode', 'i0 3E2full cathode');

set(gcf,"Position",[100 100 800 1200])

%Ds
figure(2)
t = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

nexttile
plot(soc, log10(Ds_sum(1,:)),'bo-','LineWidth',1,'MarkerFaceColor','b'); %Sum_anode
    set(gca,'Box','on',...    
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    title('3E_Simul anode')
    xlabel('soc')
    ylabel('log[Ds]')
    legend('Ds_Simul anode');

nexttile
plot(soc, log10(Ds_sum(2,:)),'bo-','LineWidth',1,'MarkerFaceColor','b'); %Sum_cathode
set(gca,'Box','on',...  
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
     title('3E_Simul cathode')
    xlabel('soc')
    ylabel('log[Ds]')
    legend('Ds_Simul cathode');

nexttile
plot(soc, log(Ds_3E2full(1,:)),'ro-','LineWidth',1,'MarkerFaceColor','r'); %3E2full_anode
set(gca,'Box','on',... 
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    hold off
    title('3_Sum anode')
    xlabel('soc')
    ylabel('log[Ds]')
    legend('Ds_Sum anode');

nexttile
plot(soc, log10(Ds_3E2full(2,:)),'ro-','LineWidth',1,'MarkerFaceColor','r'); %3E2full_cathode
set(gca,'Box','on',...   
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
   title('3E_Sum cathode')
    xlabel('soc')
    ylabel('log[Ds]')
    legend('Ds_Sum cathode');

nexttile
plot(soc, log10(Ds_sum(1,:)),'bo-','LineWidth',1,'MarkerFaceColor','b'); hold on %Sum + 3E2full anode in one fig
plot(soc, log10(Ds_3E2full(1,:)),'ro-','LineWidth',1,'MarkerFaceColor','r'); %3E2full_anode
    set(gca,'Box','on',...    
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    hold off
    title('3E_Simul + 3E_Sum Ds anode')
    xlabel('soc')
    ylabel('log[Ds]')
    legend('Ds_Simul anode', 'Ds_Sum anode');

    nexttile
plot(soc, log10(Ds_sum(2,:)),'bo-','LineWidth',1,'MarkerFaceColor','b'); hold on%Sum + 3E2full cathode in one fig
plot(soc, log10(Ds_3E2full(2,:)),'ro-','LineWidth',1,'MarkerFaceColor','r');
set(gca,'Box','on',...  
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    hold off
     title('3E_Simul + 3E_Sum Ds cathode')
    xlabel('soc')
    ylabel('log[Ds]')
    legend('Ds_Simul cathode', 'Ds_Sum cathode');

set(gcf,"Position",[920 100 800 1200])

%K_el
figure(3)
t = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile
plot(soc, K_el_sum(1,:),'bo-','LineWidth',1,'MarkerFaceColor','b'); %Sum_anode
    set(gca,'Box','on',...    
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    title('3E Simul anode')
    xlabel('soc')
    ylabel('K el')
    legend('K el Simul');


nexttile
plot(soc, K_el_3E2full(1,:),'ro-','LineWidth',1,'MarkerFaceColor','r'); %3E2full_anode
set(gca,'Box','on',... 
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    hold off
    title('3E Sum anode')
    xlabel('soc')
    ylabel('K el')
    legend('K el 3E Sum');

nexttile
plot(soc, K_el_sum(1,:),'bo-','LineWidth',1,'MarkerFaceColor','b'); hold on %Sum + 3E2full anode in one fig
plot(soc, K_el_3E2full(1,:),'ro-','LineWidth',1,'MarkerFaceColor','r'); %3E2full_anode
    set(gca,'Box','on',...    
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    hold off
    title('3E Simul + 3E Sum K el')
    xlabel('soc')
    ylabel('K el')
    legend('K el Simul', 'K el 3E Sum');

set(gcf,"Position",[100 100 400 1200])

%D_el
figure(4)
t = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile
plot(soc, log10(D_el_sum(1,:)),'bo-','LineWidth',1,'MarkerFaceColor','b'); %Sum_anode
    set(gca,'Box','on',...    
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    title('3E Simul anode')
    xlabel('soc')
    ylabel('log[D el]')
    legend('D el Simul');


nexttile
plot(soc, log10(D_el_3E2full(1,:)),'ro-','LineWidth',1,'MarkerFaceColor','r'); %3E2full_anode
set(gca,'Box','on',... 
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    hold off
    title('3E Sum anode')
    xlabel('soc')
    ylabel('log[D el]')
    legend('D el Sum');

nexttile
plot(soc, log10(D_el_sum(1,:)),'bo-','LineWidth',1,'MarkerFaceColor','b'); hold on %Sum + 3E2full anode in one fig
plot(soc, log10(D_el_3E2full(1,:)),'ro-','LineWidth',1,'MarkerFaceColor','r'); %3E2full_anode
    set(gca,'Box','on',...    
    'PlotBoxAspectRatio',[1 1 1],... 
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    grid on;
    hold off
    title('3E Simul + 3E Sum D el')
    xlabel('soc')
    ylabel('log[D el]')
    legend('D el Simul', 'D el 3 Sum');

set(gcf,"Position",[520 100 400 1200])