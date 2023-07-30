%% 1.2.1 简单的合成音乐

close all;
clear;
clc;

%% (1)

%% 生成幅度为 1 、抽样频率为 8kHz 的正弦信号

% 取一节的时间长度为 1 秒
sample_rate = 8000;
section_duration = 1;

% 产生长为 1 节 (2 拍) 的时间序列
t_double = (0 : 1/sample_rate : section_duration-1/sample_rate)';

% 1 节时间内的样点数
count_double = length(t_double);

% 确定 3 组 21 个频率值, 包含低音、中音、高音 1 2 3 4 5 6 7
Freq_0 = 349.23;
freq = Freq_0 * 2.^([0,2,4,5,7,9,11]/12);
freq_flat = freq / 2;
freq_sharp = 2 * freq;

% 产生正弦单音
double_tone_mat = sin(t_double * 2*pi*[freq_flat,freq,freq_sharp]); % 2   拍
       tone_mat = double_tone_mat(1:count_double/2,:);              % 1   拍
  half_tone_mat = double_tone_mat(1:count_double/4,:);              % 1/2 拍

tones = mat2cell([half_tone_mat;tone_mat;double_tone_mat], ...
    [count_double/4,count_double/2,count_double], ones(21,1))';

% 采用变调和时长调整标记，用于在 tones 索引
s =  -7; % sharp
f =   7; % flat
d = -21; % double
h =  21; % half

%% 用 sound 函数播放每个乐音

% 播放 6. 7. 1 2 3 4 5 6
clips = generate_clips(tones, [6-f, 7-f, 1, 2, 3, 4, 5, 6], sample_rate, 0);
sound(clips, sample_rate);

%% 拼出《东方红》片断

clips = generate_clips(tones, [5, 5-h, 6-h, 2-d, 1, 1-h, 6-f-h, 2-d], sample_rate, 0);
sound(clips, sample_rate);

%% (2)

%% 包络修正

% 迭接部分时长
overlapping_length = 0.02;

% 对 2 拍, 1 拍, 1/2 拍单音分别调整包络修正参数
double_envelope = music_envelope( ...
    [0.01,0.05,0.20,0.70,1.02], 3, [100,-100,30,-10], ...
    sample_rate);
envelope        = music_envelope( ...
    [0.01,0.05,0.15,0.30,0.52], 3, [100,-100,40,-20], ...
    sample_rate);
half_envelope   = music_envelope( ...
    [0.01,0.05,0.10,0.15,0.27], 3, [100,-100,50,-30], ...
    sample_rate);

% 绘制包络形状
close all;
subplot(3,1,1);
plot(double_envelope);
subplot(3,1,2);
plot(envelope), xlim([0,sample_rate*section_duration]);
subplot(3,1,3);
plot(half_envelope), xlim([0,sample_rate*section_duration]);

% 包含迭接部分时长，产生单音
t_double = (0 : 1/sample_rate : overlapping_length+section_duration)';
t_double = t_double(1:end-1);

count_double = length(t_double);
count        = length(0:1/sample_rate:overlapping_length+section_duration/2)-1;
count_half   = length(0:1/sample_rate:overlapping_length+section_duration/4)-1;

double_tone_mat = sin(t_double * 2*pi*[freq_flat,freq,freq_sharp]);
       tone_mat = double_tone_mat(1:count,:);
  half_tone_mat = double_tone_mat(1:count_half,:);

% 将包络与单音相乘
double_tone_mat = double_envelope .* double_tone_mat;
       tone_mat =        envelope .*        tone_mat;
  half_tone_mat =   half_envelope .*   half_tone_mat;

tones = mat2cell([half_tone_mat;tone_mat;double_tone_mat], ...
    [count_half,count,count_double],ones(21,1))';

% 播放《东方红》片段
clips = generate_clips(tones,[5,5-h,6-h,2-d,1,1-h,6-f-h,2-d], ...
    sample_rate, overlapping_length);
sound(clips,sample_rate);

%% (3)

%% 将 (2) 中的音乐升高一个八度
sound(clips, sample_rate*2);

%% 将 (2) 中的音乐降低一个八度
sound(clips, sample_rate/2);

%% 使用 resample 函数升高一个半音

clips = resample(clips,round(sample_rate*2^(-1/12)), sample_rate);
sound(clips, sample_rate);

%% (4)

% 谐波相对幅度取 0.2, 0.3
double_tone_mat = sin(t_double * 2*pi*[freq_flat,freq,freq_sharp]) ...
          + 0.2 * sin(t_double * 2*pi*[freq_flat,freq,freq_sharp]*2) ...
          + 0.3 * sin(t_double * 2*pi*[freq_flat,freq,freq_sharp]*3);
       tone_mat = double_tone_mat(1:count,:);
  half_tone_mat = double_tone_mat(1:count_half,:);

double_tone_mat = double_envelope .* double_tone_mat;
       tone_mat =        envelope .*        tone_mat;
  half_tone_mat =   half_envelope .*   half_tone_mat;

tones = mat2cell([half_tone_mat;tone_mat;double_tone_mat], ...
    [count_half,count,count_double], ones(21,1))';

% 播放《东方红》片段
clips = generate_clips(tones,[5,5-h,6-h,2-d,1,1-h,6-f-h,2-d], ...
    sample_rate, overlapping_length);
sound(clips,sample_rate);

%% (5)

clear;
clc;

%% 《东方红》

The_East_Is_Red = music_synthesis( ...
        sample_rate = 8000, ...
        speed = 'Moderato', ...
        beat_tempo = '2/4', ...
        base_frequency = 'f1', ...
        envelope_description = music_const.envelope_piano, ...
        harmonic =  music_const.harmonic_piano ...
    );

The_East_Is_Red.clips_description = [
    '5 _5 _6 ' ...
    '2 - ' ...
    '1 _1 _b6 ' ...
    '2 -'
    ];

audiowrite('The_East_Is_Red_1.wav', The_East_Is_Red.clips, The_East_Is_Red.sample_rate);

%% 《夜明》

Dawn = music_synthesis( ...
        sample_rate = 8000, ...
        speed = 'Andantino', ...
        beat_tempo = '4/4' ...
    );

clips_description_1 = [
    '0 ' ...
    '6 . _#1 7 _6 _5 ' ...
    '_5 _6 3 - _3 _5 ' ...
    '6 . _#1 7 _6 _7 ' ...
    '7 - - _3 _5 ' ...
    '6 . _#1 7 _6 _5 ' ...
    '_5 _6 3 - _2 _1 ' ...
    '_6 3 _5 _3 _2 _1 _2 ' ...
    'b6 - - - '
    ];

clips_description_2 = [
    '0 0 0 _3 _5 ' ...
    '6 . _#1 7 _6 _5 ' ...
    '_5 _6 3 - _3 _5 ' ...
    '6 . _#1 7 _6 _7 ' ...
    '7 - - _3 _5 ' ...
    '6 . _#1 7 _6 _5 ' ...
    '_5 _6 3 - _2 _1 ' ...
    '_b6 3 _5 _3 _2 _1 _2 ' ...
    'b6 - - - ' ...
    '1 - 2 _1 _2 ' ...
    '_5 __6 __5 _3 _5 3 - ' ...
    '1 - 2 _1 _2 ' ...
    '5 _6 __5 __6 3 _3 _5 ' ...
    '6 . _#1 7 _6 _5 ' ...
    '_5 _6 3 - _3 _5 ' ...
    '6 . _#1 7 _6 _7 ' ...
    '7 - - _3 _5 ' ...
    '6 . _#1 7 _6 _5 ' ...
    '_5 _6 3 - _2 _1 ' ...
    '_b6 3 _5 _3 _2 _1 _2 ' ...
    'b6 - - - '
    ];

clips_description_3 = [
    '0 0 0 _b6 _b7 ' ...
    '1 3 b7 3 ' ...
    'b6 - - _b6 _b7 ' ...
    '1 _b7 _b6 _b5 _b6 _5 _2 ' ...
    '3 - - _b6 _b7 ' ...
    '1 3 b7 3 ' ...
    'b6 b7 1 2 ' ...
    '3 - - _5 _6 ' ...
    '3 - - 2 ' ...
    '3 - 6 5 ' ...
    '3 - - - ' ...
    'b6 - - b7 ' ...
    '1 2 b7 3 ' ...
    'b6 1 b7 3 ' ...
    'b6 - - - ' ...
    'b6 - - - '
    ];

clips_description_4_1 = [
    '0 0 0 _3 _5 ' ...
    '6 . _#1 7 _6 _5 ' ...
    '_5 _6 3 - _3 _5 ' ...
    '6 . _#1 7 _6 _7 ' ...
    '7 - - _3 _5 ' ...
    '6 . _#1 7 _6 _5 ' ...
    '_5 _6 3 - _2 _1 ' ...
    '_6 3 _5 _3 _2 _1 _2 ' ...
    'b6 - - _3 _5 ' ...
    '1 - 2 _1 _2 ' ...
    '_5 __6 __5 _3 _5 3 - ' ...
    '1 - 2 _1 _2 ' ...
    '5 _6 __5 __6 3 - ' ...
    ];

clips_description_4_2 = [
    '0 0 0 3 ' ...
    '6 - - - ' ...
    '5 3 - - ' ...
    '6 - - - ' ...
    '7 - - - ' ...
    '6 - - - ' ...
    '5 3 - - ' ...
    '6 3 - - ' ...
    'b6 - - - ' ...
    '1 - - - ' ...
    '5 - - - ' ...
    '1 - - - ' ...
    '5 - - - ' ...
    ];

clips_description_5 = [
    '0 0 0 _3 _5 ' ...
    '6 . _#1 7 _6 _5 ' ...
    '_5 _6 3 - _3 _5 ' ...
    '6 . _#1 7 _6 _7 ' ...
    '7 - - _3 _5 ' ...
    '6 . _#1 7 _6 _5 ' ...
    '_5 _6 3 - _2 _1 ' ...
    '_b6 3 _5 _3 _2 _1 _2 ' ...
    'b6 - - - ' ...
    'b6 - - _2 _1 ' ...
    '_b6 3 _5 _3 _2 _1 _2 ' ...
    'b6 - - - ' ...
    ];

Dawn.base_frequency = 'bb';
Dawn.envelope_description = music_const.envelope_piano;
Dawn.harmonic = music_const.harmonic_piano;
Dawn.clips_description = clips_description_1;
clips_1 = Dawn.clips;

Dawn.base_frequency = 'bb1';
Dawn.envelope_description = music_const.envelope_shakuhachi;
Dawn.harmonic = music_const.harmonic_shakuhachi;
Dawn.clips_description = clips_description_2;
clips_2 = Dawn.clips;

Dawn.base_frequency = 'bb';
Dawn.envelope_description = music_const.envelope_piano;
Dawn.harmonic = music_const.harmonic_piano;
Dawn.clips_description = clips_description_3;
clips_3 = Dawn.clips;

Dawn.base_frequency = 'bb1';
Dawn.envelope_description = music_const.envelope_shakuhachi;
Dawn.harmonic = music_const.harmonic_shakuhachi;
Dawn.clips_description = clips_description_4_1;
clips_4_1 = Dawn.clips;

Dawn.base_frequency = 'bb';
Dawn.envelope_description = music_const.envelope_piano;
Dawn.harmonic = music_const.harmonic_piano;
Dawn.clips_description = clips_description_4_2;
clips_4_2 = Dawn.clips;
length_difference = length(clips_4_1) - length(clips_4_2);
if length_difference < 0
    clips_4 = 0.6*[clips_4_1; zeros(-length_difference, 1)] + 0.4*clips_4_2;
else
    clips_4 = 0.6*clips_4_1 + 0.4*[clips_4_2; zeros(length_difference, 1)];
end

Dawn.clips_description = clips_description_5;
clips_5 = Dawn.clips;

audiowrite('Dawn.wav', [clips_1; clips_2; clips_3; clips_4; clips_5], Dawn.sample_rate);

%% 《玛丽有只小羊羔》

Mary_Had_a_Little_Lamb = music_synthesis( ...
        sample_rate = 8000, ...
        speed = 96, ...
        beat_tempo = '4/4', ...
        base_frequency = 'c1', ...
        envelope_description = music_const.envelope_piano, ...
        harmonic = music_const.harmonic_piano ...
    );

Mary_Had_a_Little_Lamb.clips_description = [
    '3 2 1 2 ' ...
    '3 3 3 - ' ...
    '2 2 2 - ' ...
    '3 3 3 - ' ...
    '3 2 1 2 ' ...
    '3 3 3 1 ' ...
    '2 2 3 2 ' ...
    '1 - - - ' ...
];

audiowrite('Mary_Had_a_Little_Lamb.wav', ...
    Mary_Had_a_Little_Lamb.clips, Mary_Had_a_Little_Lamb.sample_rate);

%% 1.2.2 用傅里叶级数分析音乐

clear;
clc;

%% 载入 guitar.mat

load('Guitar.MAT');

%% 用 wavread 函数载入 fmt.wav 文件

[clips, sample_rate] = audioread('fmt.wav');
sound(clips, sample_rate);

%% (7)

%% 去除非线性谐波和噪声

sample_rate = 8000;

% 估计基音周期
R_realwave = xcorr(realwave);                  % 自相关函数
[peak_y, peak_x] = findpeaks(R_realwave);      % 包含全部峰值
peak_x = peak_x (peak_y > 0.75 * max(peak_y)); % 仅保留主要峰值

point_count = floor(length(peak_x)/2);         % 逐差点数
point_distance = ceil(length(peak_x)/2);       % 逐差点间距
approx_T = (sum(peak_x(end-point_count+1:end)) - sum(peak_x(1:point_count))) ...
    / point_distance / point_count;

% 确定基音周期
resample_wave = resample(realwave, 100, 1);
R_resample_wave = xcorr(resample_wave);
[peak_y, peak_x] = findpeaks(R_resample_wave, MinPeakDistance=approx_T*100*0.8);

point_count = floor(length(peak_x)/2);
point_distance = ceil(length(peak_x)/2);
resample_wave_T = (sum(peak_x(end-point_count+1:end)) - sum(peak_x(1:point_count))) ...
    / point_distance / point_count;            % 重采样波形一周期样点数

% 时间同步平均
resample_wave_length = length(resample_wave);                          % 重采样波形样点数
resample_wave_T_count = floor(resample_wave_length / resample_wave_T); % 重采样波形周期数
resample_wave_single_T = tsa(resample_wave, 1, 0:resample_wave_T:resample_wave_length);

% 周期重复
resample_wave2proc_2 = repmat(resample_wave_single_T, resample_wave_T_count, 1);
resample_wave2proc_2 = [resample_wave2proc_2; ...
    resample_wave_single_T(1 : resample_wave_length-length(resample_wave2proc_2))];

% 还原采样率
wave2proc_2 = resample(resample_wave2proc_2, 1, 100);

%% 绘制波形, 自相关函数, 峰值

close all;
subplot(4,1,1);
plot(realwave);
title('原波形');
xlim([-243,243]);

subplot(4,1,2);
plot(resample_wave);
title('重采样波形');
xlim([-24300,24300]);

subplot(4,1,3);
plot(R_resample_wave);
title('自相关函数');
xlim([0,48599]);

subplot(4,1,4);
plot(peak_x, peak_y,'x');
title('自相关函数峰值');
xlim([0,48599]);

%% 绘制 wave2proc

close all;
subplot(2,1,1);
plot(wave2proc);
title('wave2proc');
subplot(2,1,2);
plot(wave2proc_2);
title('wave2proc 2');

%% (8)

%% 基频

base_frequency = 1 / (resample_wave_T / 100 / sample_rate);

%% 用傅里叶级数 DFS 分析谐波分量

% 近似取一个周期求傅里叶级数
wave2proc_T = round(resample_wave_T / 100);       % wave2proc 一周期样点数
wave2proc_single_T = wave2proc(1 : wave2proc_T);  % wave2proc 单周期
x_shift = (-wave2proc_T/2 : wave2proc_T/2-1)';    % 零点居中 x 坐标
DFS = abs(fftshift(fft(wave2proc_single_T))) / length(wave2proc_single_T);
stem(x_shift, DFS / DFS(x_shift == 1)); % 基频幅值调整为 1
title('傅里叶级数');

%% 用傅里叶变换 DFT 分析谐波分量

wave2proc_length = length(wave2proc);
x_shift = sample_rate / wave2proc_length * ...
    (-floor(wave2proc_length/2) : ceil(wave2proc_length/2-1))';
DFT = abs(fftshift(fft(wave2proc)));
base_frequency_peak = max(DFT(x_shift-base_frequency<10 & x_shift-base_frequency>-10));
plot(x_shift, DFT / base_frequency_peak);
title('傅里叶变换');

%% 提高分辨率

x_shift = sample_rate * 100 / resample_wave_length * ...
    (-floor(resample_wave_length/2) : ceil(resample_wave_length/2-1))';
DFT = abs(fftshift(fft(resample_wave2proc_2)));
base_frequency_peak = max(DFT(x_shift-base_frequency<10 & x_shift-base_frequency>-10));
plot(x_shift, DFT / base_frequency_peak);
xlim([-4000, 4000]);
title('提高分辨率傅里叶变换');

%% 时域重复 20 周期后傅里叶变换

repeat_resample_wave = repmat(resample_wave_single_T, 20, 1);
repeat_resample_wave_length = length(repeat_resample_wave);
x_shift = sample_rate * 100 / repeat_resample_wave_length * ...
    (-floor(repeat_resample_wave_length/2) : ceil(repeat_resample_wave_length/2-1))';
DFT = abs(fftshift(fft(repeat_resample_wave)));
base_frequency_peak = max(DFT(x_shift-base_frequency<10 & x_shift-base_frequency>-10));
plot(x_shift, DFT / base_frequency_peak);
xlim([-4000, 4000]);
title('时域重复 20 周期后傅里叶变换');

%% 时域重复 100 周期后傅里叶变换

repeat_resample_wave = repmat(resample_wave_single_T, 100, 1);
repeat_resample_wave_length = length(repeat_resample_wave);
x_shift = sample_rate * 100 / repeat_resample_wave_length * ...
    (-floor(repeat_resample_wave_length/2) : ceil(repeat_resample_wave_length/2-1))';
DFT = abs(fftshift(fft(repeat_resample_wave)));
base_frequency_peak = max(DFT(x_shift-base_frequency<10 & x_shift-base_frequency>-10));
plot(x_shift, DFT / base_frequency_peak);
xlim([-4000, 4000]);
title('时域重复 100 周期后傅里叶变换');

%% (9)

close all;
clear;
clc;

[guitar_clips, sample_rate] = audioread('fmt.wav');

guitar_analysis = music_analysis( ...
    clips = guitar_clips ...
    );

guitar_analysis.analyse();
disp(guitar_analysis.tone);
disp(guitar_analysis.beat);

guitar_synthesis = music_synthesis( ...
        sample_rate = 8000, ...
        speed = 'Moderato', ...
        beat_tempo = '4/4', ...
        base_frequency = 'c1', ...
        envelope_description = music_const.envelope_guitar, ...
        harmonic = music_const.harmonic_guitar ...
    );

guitar_synthesis.clips_description = [
    'b5 - - - ' ...
    'b7 b6 2 3 b5 b6 b4 ' ...
    '2 - __3 __b6 __b7 ' ...
    '3 - 3 - ' ...
    'b6 - _3 _b6 6 ' ...
    'b6 _5 _4 _3 _2 1 ' ...
    'b7 __7 __b7 2 _1 ' ...
    '__b7 _1 b7 b6 _b7 __2 ' ...
    '_b7 __b5 b6 - b5 _b6 1 - ' ...
    ];

audiowrite('guitar_synthesis.wav', guitar_synthesis.clips, guitar_synthesis.sample_rate);

%% (10)

%% 《东方红》

The_East_Is_Red = music_synthesis( ...
        sample_rate = 8000, ...
        speed = 'Moderato', ...
        beat_tempo = '2/4', ...
        base_frequency = 'f1', ...
        envelope_description = music_const.envelope_guitar, ...
        harmonic =  music_const.harmonic_guitar ...
    );

The_East_Is_Red.clips_description = [
    '5 _5 _6 ' ...
    '2 - ' ...
    '1 _1 _b6 ' ...
    '2 -'
    ];

audiowrite('The_East_Is_Red_2.wav', The_East_Is_Red.clips, The_East_Is_Red.sample_rate);

%% (11)

% C 调谐波分量
harmonic_guitar_c = guitar_analysis.harmonic.c1;

% D 调谐波分量
harmonic_guitar_d = guitar_analysis.harmonic.d1;

% F 调谐波分量
harmonic_guitar_f = guitar_analysis.harmonic.f1;

% G 调谐波分量
harmonic_guitar_g = guitar_analysis.harmonic.g1;

The_East_Is_Red = music_synthesis( ...
        sample_rate = 8000, ...
        speed = 'Moderato', ...
        beat_tempo = '2/4', ...
        base_frequency = 'f1', ...
        envelope_description = music_const.envelope_guitar ...
    );
The_East_Is_Red.harmonic = harmonic_guitar_c;
The_East_Is_Red.clips_description = '5 _5';
clips_1 = The_East_Is_Red.clips;

The_East_Is_Red.harmonic = harmonic_guitar_d;
The_East_Is_Red.clips_description = '_6';
clips_2 = The_East_Is_Red.clips;

The_East_Is_Red.harmonic = harmonic_guitar_g;
The_East_Is_Red.clips_description = '2 -';
clips_3 = The_East_Is_Red.clips;

The_East_Is_Red.harmonic = harmonic_guitar_f;
The_East_Is_Red.clips_description = '1 _1';
clips_4 = The_East_Is_Red.clips;

The_East_Is_Red.harmonic = harmonic_guitar_d;
The_East_Is_Red.clips_description = '_b6';
clips_5 = The_East_Is_Red.clips;

clips = [clips_1; clips_2; clips_3; clips_4; clips_5; clips_3];

audiowrite('The_East_Is_Red_3.wav', clips, The_East_Is_Red.sample_rate);
