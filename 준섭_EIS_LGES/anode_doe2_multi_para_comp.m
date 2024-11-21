clc, clear, close all;

% load paras_integ_dists respectively first.
addpath('G:\공유 드라이브\Battery Software Lab\Projects\LGES 2023\발표 및 공유자료\10월 미팅자료\흑연 2종 결과 데이터\multi-soc_3E_Simul\blending')
addpath('G:\공유 드라이브\Battery Software Lab\Projects\LGES 2023\발표 및 공유자료\10월 미팅자료\흑연 2종 결과 데이터\multi-soc_3E_Simul\천연_100')
addpath('G:\공유 드라이브\Battery Software Lab\Projects\LGES 2023\발표 및 공유자료\10월 미팅자료\기존_3E_SOC 통합 데이터\3E_Simul')



para_integ_dist = load('3E simul P2D + soc integrated + DRT + DDT para_sum_integ_dist.mat');
para_blend = load('인조_천연 Blending음극_#1 3E simul P2D + soc integrated + DRT + DDT para_sum_integ_dist.mat');
para_nat = load('천역흑연100_음극_#2 3E simul P2D + soc integrated + DRT + DDT para_sum_integ_dist.mat');


soc = 20:10:80;
%multi data
multi_i0a = para_integ_dist.paras_integ_dist(3,:);
multi_Cdla = para_integ_dist.paras_integ_dist(4,:);
multi_Dsa = para_integ_dist.paras_integ_dist(5,:);
multi_Ava = para_integ_dist.paras_integ_dist(6,:);

multi_i0p = para_integ_dist.paras_integ_dist(8,:);
multi_Cdlp = para_integ_dist.paras_integ_dist(9,:);
multi_Dsp = para_integ_dist.paras_integ_dist(10,:);
multi_Avp = para_integ_dist.paras_integ_dist(11,:);

multi_Kel = para_integ_dist.paras_integ_dist(16,:);
multi_Del = para_integ_dist.paras_integ_dist(17,:);

multi_Std_DRTa = para_integ_dist.paras_integ_dist(12,:);
multi_Std_DRTc = para_integ_dist.paras_integ_dist(14,:);

multi_data = [multi_i0a; multi_Cdla; multi_Dsa; multi_Ava; multi_i0p;multi_Cdlp;multi_Dsp;multi_Avp;
    multi_Kel;multi_Del;multi_Std_DRTa;multi_Std_DRTc];

% blend data
i0a_blend = para_blend.paras_integ_dist(3,:);
Cdla_blend = para_blend.paras_integ_dist(4,:);
Dsa_blend = para_blend.paras_integ_dist(5,:);
Ava_blend = para_blend.paras_integ_dist(6,:);

i0p_blend = para_blend.paras_integ_dist(8,:);
Cdlp_blend = para_blend.paras_integ_dist(9,:);
Dsp_blend = para_blend.paras_integ_dist(10,:);
Avp_blend = para_blend.paras_integ_dist(11,:);

Kel_blend = para_blend.paras_integ_dist(16,:);
Del_blend = para_blend.paras_integ_dist(17,:);

Std_DRTa_blend = para_blend.paras_integ_dist(12,:);
Std_DRTc_blend = para_blend.paras_integ_dist(14,:);

blend_data = [i0a_blend; Cdla_blend; Dsa_blend; Ava_blend; i0p_blend; Cdlp_blend; Dsp_blend; Avp_blend; Kel_blend; Del_blend;Std_DRTa_blend; Std_DRTc_blend ];

% natural data
i0a_nat = para_nat.paras_integ_dist(3,:);
Cdla_nat = para_nat.paras_integ_dist(4,:);
Dsa_nat = para_nat.paras_integ_dist(5,:);
Ava_nat = para_nat.paras_integ_dist(6,:);

i0p_nat = para_nat.paras_integ_dist(8,:);
Cdlp_nat = para_nat.paras_integ_dist(9,:);
Dsp_nat = para_nat.paras_integ_dist(10,:);
Avp_nat = para_nat.paras_integ_dist(11,:);

Kel_nat = para_nat.paras_integ_dist(16,:);
Del_nat = para_nat.paras_integ_dist(17,:);

Std_DRTa_nat = para_nat.paras_integ_dist(12,:);
Std_DRTc_nat = para_nat.paras_integ_dist(14,:);

natural_data = [i0a_nat; Cdla_nat; Dsa_nat; Ava_nat; i0p_nat; Cdlp_nat; Dsp_nat; Avp_nat; Kel_nat; Del_nat;Std_DRTa_nat;Std_DRTc_nat];

% plot to compare each parameter 
name = {'i0a' 'Cdla' 'Dsa' 'Ava' 'i0p' 'Cdlp' 'Dsp' 'Avp' 'Kel' 'Del' 'Std DRTa' 'Std DRTc'};

c_mat_single = lines(length(name));
c_mat_multi = jet(length(name));

for i = 1:length(name)
    figure
    hold on
    plot(soc, multi_data(i,:),'Color',c_mat_multi(9,:),'LineStyle','--','Marker','o','MarkerFaceColor',c_mat_multi(9,:));
    plot(soc, blend_data(i,:),'Color',c_mat_single(8,:),'LineStyle','--','Marker','o','MarkerFaceColor',c_mat_single(8,:)); hold on;
    plot(soc, natural_data(i,:),'Color',c_mat_single(2,:),'LineStyle','--','Marker','o','MarkerFaceColor',c_mat_single(2,:));

    hold off
    title(['parameter ' name{i} ' comparison'])
    legend(['Multi-soc ' name{i}],['Blend ' name{i}],['Natural ' name{i}])
    box on;
end

