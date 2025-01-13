clc, clear, close all;

Refs = {
    [57 207 134 150 195],
    [136 115 149 190 57],
    [50 219 203 233],
    [257 228 179 175],
    [235 59 148 153],
    [124 80 89 140],
    [19 119 143 161],
    [269 50 141 252]
};

figure()
t = tiledlayout(3,3,"TileSpacing","compact","Padding","compact");

for i = 1:length(Refs)
    nexttile
    plot(Refs{i},'o','MarkerFaceColor','black','MarkerSize',8) 

    set(gca,'PlotBoxAspectRatio',[1 1 1], 'fontunit', 'point', 'fontsize', 11, 'fontname', 'TimesNewRoman')
    ylim([0 300])
    grid on;
    box on;
    xlabel('Sample #')
    ylabel('Cycle')
    title(sprintf('Ref%d', i)) 
end

set(gcf,'Position',[100 100 1200 1200])