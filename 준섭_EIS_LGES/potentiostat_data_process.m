close all; clear; clc;

file_folder = 'G:\공유 드라이브\BSL-Data\카이스트_단락셀\카이스트 단락셀\2차 셀 데이터\RPT\RPT';
file_name = 'Ref_RPT(o)_1_01_DC.txt'
data_folder = fullfile(file_folder, file_name);

data = readtable(data_folder, 'NumHeaderLines', 23, ReadVariableNames=false);

% time_str을 테이블에서 추출

time_str = data.Var2; % 시간 문자열 열 추출
I = data.Var9;
V = data.Var10;

% 형식 갯수 파악 
lengths = cellfun(@length, time_str);

[unique_lengths, ~, idx] = unique(lengths);
counts = histc(idx, 1:numel(unique_lengths));

for i = 1:length(unique_lengths)
    fprintf('길이 %d인 샘플의 개수: %d\n\n\n\n', unique_lengths(i), counts(i));
end

% hh:mm:ss 형식과 mm:ss 형식을 구분
hhmmssSSS_idx = cellfun(@(x) length(x) == 12, time_str);
hhmmssSS_idx = cellfun(@(x) length(x) == 11, time_str);
hhmmssS_idx = cellfun(@(x) length(x) == 10, time_str);
mmssSSS_idx = cellfun(@(x) length(x) == 9, time_str);
mmss_idx = cellfun(@(x) length(x) == 5, time_str);
mss_idx = cellfun(@(x) length(x) == 4, time_str);
ss_idx = cellfun(@(x) length(x) == 2, time_str);

% string 나누기
hhmmssSSS_str = time_str(hhmmssSSS_idx); %수정 불필요
hhmmssSS_str = time_str(hhmmssSS_idx); %0 맨 뒤 추가
hhmmssS_str = time_str(hhmmssS_idx); % 00맨 뒤 추가
mmssSSS_str = time_str(mmssSSS_idx); % 00맨 앞 추가
mmss_str = time_str(mmss_idx); % 00맨 앞, 000맨 뒤 추가
mss_str = time_str(mss_idx); %00 0맨 앞, 000맨 뒤 추가
ss_str = time_str(ss_idx); %00 00맨 앞, 000맨 뒤 추가

% hh:mm:ss.SSS 형식 맞추기 위해 0 추가된 셀 생성

m_hhmmssSSS_str = cellfun(@(x) [x], hhmmssSSS_str, 'UniformOutput', false);
m_hhmmssSS_str = cellfun(@(x) [x '0'], hhmmssSS_str, 'UniformOutput', false);
m_hhmmssS_str = cellfun(@(x) [x '00'], hhmmssS_str, 'UniformOutput', false);
m_mmssSSS_str = cellfun(@(x) ['00:' x], mmssSSS_str, 'UniformOutput', false);
m_mmss_str = cellfun(@(x) ['00:' x '.000'], mmss_str, 'UniformOutput', false);
m_mss_str = cellfun(@(x) ['00:0' x '.000'], mss_str, 'UniformOutput', false);
m_ss_str = cellfun(@(x) ['00:00:' x '.000'], ss_str, 'UniformOutput', false);

% Duration 변환
hhmmssSSS_duration = duration(m_hhmmssSSS_str,'InputFormat','hh:mm:ss.SSS');
hhmmssSS_duration = duration(m_hhmmssSS_str,'InputFormat','hh:mm:ss.SSS');
hhmmssS_duration = duration(m_hhmmssS_str,'InputFormat','hh:mm:ss.SSS');
mmssSSS_duration = duration(m_mmssSSS_str,'InputFormat','hh:mm:ss.SSS');
mmss_duration = duration(m_mmss_str,'InputFormat','hh:mm:ss.SSS');
mss_duration = duration(m_mss_str,'InputFormat','hh:mm:ss.SSS');
ss_duration = duration(m_ss_str,'InputFormat','hh:mm:ss.SSS');

% Duration 합치기
ini_zero = zeros(length(time_str),3);
time = duration(ini_zero); %초기화된 duration 배열 생성

time(hhmmssSSS_idx) = hhmmssSSS_duration;
time(hhmmssSS_idx) = hhmmssSS_duration;
time(hhmmssS_idx) = hhmmssS_duration;
time(mmssSSS_idx) = mmssSSS_duration;
time(mmss_idx) = mmss_duration;
time(mss_idx) = mss_duration;
time(ss_idx) = ss_duration;

%plot
figure()
hold on

yyaxis left
plot(time,I,'b-', 'LineWidth', 1)
ylabel('Current [A]')


yyaxis right
plot(time,V,'r-', 'LineWidth', 1)
ylabel('Voltage [V]')
hold off

set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')

title ('I/V vs time graph')
xlabel('Time [s]')

legend ('Current', 'Voltage')

%Cycle 별 확대 플랏
cycle = data.Var5;

for i = 1:max(cycle);

idx = (cycle == i);

figure()
hold on

yyaxis left
plot(time(idx), I(idx), 'b-', 'LineWidth', 1, 'DisplayName', 'Current(A)');
ylabel('Current (A)');

yyaxis right
plot(time(idx), V(idx), 'r-', 'LineWidth', 1, 'DisplayName', 'Voltage(V)');
ylabel('Voltage (V)');

xlabel('Time (s)');
title('Time vs Current/Voltage for ');
legend(['cycle', num2str(i)], 'Location','southwest');
grid on;
box on;

end