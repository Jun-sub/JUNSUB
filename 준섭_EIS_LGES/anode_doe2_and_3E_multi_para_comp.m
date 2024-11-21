clc, clear, close all;

% load para_sums respectively first.
addpath('G:\공유 드라이브\Battery Software Lab\Projects\LGES 2023\발표 및 공유자료\10월 미팅자료\흑연 2종 결과 데이터\인조_천연 혼합\Full')
addpath('G:\공유 드라이브\Battery Software Lab\Projects\LGES 2023\발표 및 공유자료\10월 미팅자료\흑연 2종 결과 데이터\천연_100%\Full')
addpath('G:\공유 드라이브\Battery Software Lab\Projects\LGES 2023\발표 및 공유자료\10월 미팅자료\기존_3E_Dist_integ 데이터\Full')



para_3E = load('full P2D + DRT + DDT para_sum.mat');
para_blend = load('인조_천연 Blending음극_#1 Full P2D + DRT + DDT para_sum.mat');
para_nat = load('천역흑연100_음극_#2 Full P2D + DRT + DDT para_sum.mat');


soc = 20:10:80;

%3E data
i0 = para_3E.para_sum(12,:);
Cdl = para_3E.para_sum(13,:);
Ds = para_3E.para_sum(14,:);
Av = para_3E.para_sum(17,:);

Kel = para_3E.para_sum(15,:);
Del = para_3E.para_sum(16,:);

Std_DRT = para_3E.para_sum(18,:);
Std_DDT = para_3E.para_sum(19,:);

data = [i0; Cdl; Ds; Av;Kel;Del;Std_DRT;Std_DDT];
data = data(:,3:9);
% blend data
i0_blend = para_blend.para_sum(12,:);
Cdl_blend = para_blend.para_sum(13,:);
Ds_blend = para_blend.para_sum(14,:);
Av_blend = para_blend.para_sum(17,:);

Kel_blend = para_blend.para_sum(15,:);
Del_blend = para_blend.para_sum(16,:);

Std_DRT_blend = para_blend.para_sum(18,:);
Std_DDT_blend = para_blend.para_sum(19,:);

blend_data = [i0_blend; Cdl_blend; Ds_blend; Av_blend; Kel_blend; Del_blend;Std_DRT_blend; Std_DDT_blend ];
blend_data = blend_data(:,3:9);
% natural data
i0_nat = para_nat.para_sum(12,:);
Cdl_nat = para_nat.para_sum(13,:);
Ds_nat = para_nat.para_sum(14,:);
Av_nat = para_nat.para_sum(17,:);

Kel_nat = para_nat.para_sum(15,:);
Del_nat = para_nat.para_sum(16,:);

Std_DRT_nat = para_nat.para_sum(18,:);
Std_DDT_nat = para_nat.para_sum(19,:);

natural_data = [i0_nat; Cdl_nat; Ds_nat; Av_nat; Kel_nat; Del_nat;Std_DRT_nat;Std_DDT_nat];
natural_data = natural_data(:,3:9);
% plot to compare each parameter 
name = {'i0' 'Cdl' 'Ds' 'Av' 'Kel' 'Del' 'Std DRT' 'Std DDT'};

c_mat_single = lines(10);
c_mat_3E = jet(10);

for i = 1:length(name)
    figure
    hold on
    plot(soc, data(i,:),'Color',c_mat_3E(9,:),'LineStyle','--','Marker','o','MarkerFaceColor',c_mat_3E(9,:));
    plot(soc, blend_data(i,:),'Color',c_mat_single(8,:),'LineStyle','--','Marker','o','MarkerFaceColor',c_mat_single(8,:)); hold on;
    plot(soc, natural_data(i,:),'Color',c_mat_single(2,:),'LineStyle','--','Marker','o','MarkerFaceColor',c_mat_single(2,:));

    hold off
    title(['parameter ' name{i} ' comparison'])
    legend(['3E ' name{i}],['Blend ' name{i}],['Natural ' name{i}])
    box on;
end

