clc, clear, close all;

% load paras_integ_dist & para_sum respectively first.
addpath('G:\공유 드라이브\Battery Software Lab\Projects\LGES 2023\발표 및 공유자료\10월 미팅자료\기존_3E_Dist_integ 데이터\3E_Simul')
addpath('G:\공유 드라이브\Battery Software Lab\Projects\LGES 2023\발표 및 공유자료\10월 미팅자료\기존_3E_SOC 통합 데이터\3E_Simul')

para_sum = load('3E simul P2D + DRT + DDT para_sum.mat');
para_integ_dist = load('3E simul P2D + soc integrated + DRT + DDT para_sum_integ_dist.mat');

soc = 20:10:80;

%single data
single_i0a = para_sum.para_sum(16,3:9);
single_Cdla = para_sum.para_sum(17,3:9);
single_Dsa = para_sum.para_sum(18,3:9);
single_Ava = para_sum.para_sum(19,3:9);

single_i0p = para_sum.para_sum(21,3:9);
single_Cdlp = para_sum.para_sum(22,3:9);
single_Dsp = para_sum.para_sum(23,3:9);
single_Avp = para_sum.para_sum(24,3:9);

single_Kel = para_sum.para_sum(25,3:9);
single_Del = para_sum.para_sum(26,3:9);
single_Del(:,end) = single_Del(:,end-1);

single_data = [single_i0a; single_Cdla; single_Dsa; single_Ava; single_i0p;single_Cdlp;single_Dsp;single_Avp;
    single_Kel;single_Del];
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

multi_data = [multi_i0a; multi_Cdla; multi_Dsa; multi_Ava; multi_i0p;multi_Cdlp;multi_Dsp;multi_Avp;
    multi_Kel;multi_Del];
% plot to compare each parameter 
name = {'i0a' 'Cdla' 'Dsa' 'Ava' 'i0p' 'Cdlp' 'Dsap' 'Avp' 'Kel' 'Del'};

c_mat_single = lines(length(name));
c_mat_multi = jet(length(name));

% multiple plot (all)
for i = 1:length(name)
figure
plot(soc, single_data(i,:),'Color',c_mat_single(9,:),'LineStyle','--','Marker','o','MarkerFaceColor',c_mat_single(9,:)); hold on;
plot(soc, multi_data(i,:),'Color',c_mat_multi(9,:),'LineStyle','--','Marker','o','MarkerFaceColor',c_mat_multi(9,:));
hold off
title(['parameter ' name{i} ' comparison'])
legend(['single ' name{i}],['multi ' name{i}])
box on;
end
