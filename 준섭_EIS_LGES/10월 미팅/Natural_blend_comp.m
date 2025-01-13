clc, clear, close all;

para_blend = load("C:\Users\admin\Documents\GitHub\JunSub\준섭_EIS_LGES\10월 미팅\Av_integ_V2\241206_Blend\3E simul P2D + soc integrated + DRT + DDT para_sum_integ_dist.mat");
para_natural = load("C:\Users\admin\Documents\GitHub\JunSub\준섭_EIS_LGES\10월 미팅\Av_integ_V2\241207_Natural\3E simul P2D + soc integrated + DRT + DDT para_sum_integ_dist.mat");

para_blend = para_blend.paras_integ_dist;
para_natural = para_natural.paras_integ_dist;

soc_vec = 20:10:80;
col_vec = [8 9 10 11 16 17];

anode_name = {'Blending' 'Natural'};
para_name = {'i0' 'C_dl' 'Ds' 'Av' 'Kel' 'Del'};

figure()
t = tiledlayout(2,3,"TileSpacing","compact","Padding","compact");

    for k = 1:length(col_vec)
        nexttile
        hold on;

        plot(soc_vec,para_blend(col_vec(k),:),'Marker','o','MarkerEdgeColor','auto','MarkerFaceColor','auto','LineWidth',1.2)
        plot(soc_vec,para_natural(col_vec(k),:),'Marker','o','MarkerEdgeColor','auto','MarkerFaceColor','auto','LineWidth',1.2)
        
        hold off;

        title(para_name(k))
        legend(anode_name,'Location','best')
        xlabel('SOC')
        grid on;
        box on;
        set(gca,'PlotBoxAspectRatio',[1 1 1],'FontUnits','points','FontSize',10,'FontName','TimesNewRoman')
    end 

    set(gcf,"Position",[300 100 1200 800])
