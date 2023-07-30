classdef music_analysis < handle
    
    properties (SetObservable, AbortSet)
        sample_rate = 8000   % 采样率
        resample_times = 100 % 重采样倍数

        STFT_window_width = 512     % 确定节拍时 STFT 海明窗宽度
        STFT_noverlap = 0.95        % 确定节拍时 STFT 重叠长度 (百分比)
        dot_step = 5                % 估计节拍时计算频率向量内积所取样点间隔
        dot_min_peak_height = 0.1   % 估计节拍时所取频率向量内积最小峰值
        dot_min_peak_distance = 0.1 % 估计节拍时所取频率向量内积峰值间距 (秒), 应小于节拍时长

        T_analyse_length = 400         % 分析频率时所取样点数
        T_DFT_min_peak = 0.1           % 估计基音周期时 允许基频 DFT 低于谐波 (百分比)
        T_time_threshold = 0.80        % 确定基音周期时 相邻峰时间阈值 (百分比)
        DFT_repeat_times = 100         % 确定谐波分量时 DFT 重复周期
        DFT_max_frequency = 4000       % 确定谐波分量时 频率范围
        DFT_base_frequency_range = 10  % 确定谐波分量时 寻找基频的范围 (1/基音周期 +- range)
        DFT_min_harmonic = 0.01        % 确定谐波分量时 所取谐波分量最小值 (百分比)

        beat_volume_threshold_next = 0.8 % 确定节奏时 认为是下一音调的音量阈值

        clips = [] % 待分析片段
    end

    properties (SetAccess = private)
        divide_time % 分隔时间
        middle_time % 取样时间
        beat        % 节拍
        tone        % 音调
        harmonic    % 谐波分量
    end
    
    methods % 初始化方法
        function obj = music_analysis(varargin)
            n = nargin;
            for i = 1:2:n
                varable_name_cell  = varargin(i);
                varable_name       = varable_name_cell{1};
                varable_cell       = varargin(i+1);
                varable            = varable_cell{1};
                obj.(varable_name) = varable;
            end
            addlistener(obj, 'sample_rate',                'PostSet', @music_analysis.set_callback);
            addlistener(obj, 'resample_times',             'PostSet', @music_analysis.set_callback);
            addlistener(obj, 'STFT_window_width',          'PostSet', @music_analysis.set_callback);
            addlistener(obj, 'STFT_noverlap',              'PostSet', @music_analysis.set_callback);
            addlistener(obj, 'dot_step',                   'PostSet', @music_analysis.set_callback);
            addlistener(obj, 'dot_min_peak_height',        'PostSet', @music_analysis.set_callback);
            addlistener(obj, 'dot_min_peak_distance',      'PostSet', @music_analysis.set_callback);
            addlistener(obj, 'T_analyse_length',           'PostSet', @music_analysis.set_callback);
            addlistener(obj, 'T_DFT_min_peak',             'PostSet', @music_analysis.set_callback);
            addlistener(obj, 'T_time_threshold',           'PostSet', @music_analysis.set_callback);
            addlistener(obj, 'DFT_repeat_times',           'PostSet', @music_analysis.set_callback);
            addlistener(obj, 'DFT_max_frequency',          'PostSet', @music_analysis.set_callback);
            addlistener(obj, 'DFT_base_frequency_range',   'PostSet', @music_analysis.set_callback);
            addlistener(obj, 'DFT_min_harmonic',           'PostSet', @music_analysis.set_callback);
            addlistener(obj, 'beat_volume_threshold_next', 'PostSet', @music_analysis.set_callback);
            addlistener(obj, 'clips',                      'PostSet', @music_analysis.set_callback);
        end
    end
        
    methods
        function analyse(obj)
            try
                obj.divide();
            catch
                error('music_analysis: 时间分隔出错');
            end
            obj.tone = cell(length(obj.middle_time), 1);
            obj.harmonic = cell(length(obj.middle_time), 1);
            for i = 1 : length(obj.middle_time)
                begin_point = round(obj.middle_time(i)*obj.sample_rate)-obj.T_analyse_length/2;
                end_point   = round(obj.middle_time(i)*obj.sample_rate)+obj.T_analyse_length/2;
                clips_i = obj.clips(max(1, begin_point) : min(end_point, end));
                try
                    [frequency_i, harmonic_i] = obj.frequency_harmonic(clips_i);
                catch
                    error(['music_analysis: 片段', num2str(i), '频率分析出错']);
                end
                obj.tone{i} = music_analysis.tone_recognition(frequency_i);
                obj.harmonic{i} = [1; harmonic_i];
            end
            try
                obj.beat_analyse();
            catch
                error('music_analysis: 节奏分析出错');
            end
        end

        function divide(obj)
            % 时间分隔

            window   = hamming(obj.STFT_window_width);
            noverlap = round(obj.STFT_window_width * obj.STFT_noverlap);
            noverlap = min(noverlap, obj.STFT_window_width-1);
            nfft     = 2^nextpow2(obj.STFT_window_width);
            
            % STFT
            [S, F, T] = spectrogram(obj.clips, window, noverlap, nfft, obj.sample_rate);
            S = abs(S);
            figure;
            subplot(4,1,1); % 绘制时频图
            mesh(T, F, mag2db(S));
            view(0, 90);
            xlim([T(1), T(end)]);
            ylabel('频率 / Hz');
            
            % 比较频率向量相似度
            frequency_dot_2 = zeros(length(T)-obj.dot_step, 1);
            for i = 1 : length(T)-obj.dot_step
                vector_1 = S(:, i);
                vector_2 = S(:, i+obj.dot_step);
                frequency_dot_2(i) = dot(vector_1, vector_2)^2 / ...
                    (dot(vector_1, vector_1) * dot(vector_2, vector_2));
            end
            subplot(4,1,2); % 绘制内积平方图
            plot(T(1:end-obj.dot_step), 1 - frequency_dot_2);
            xlim([T(1), T(end)]);
            hold on;
            
            [~, peak_x] = findpeaks(1-frequency_dot_2, ...
                MinPeakHeight = obj.dot_min_peak_height, ...
                MinPeakDistance = sum(T<obj.dot_min_peak_distance));
            obj.divide_time = T(peak_x)';                          % 峰值出现 (节奏改变) 时间
            plot(obj.divide_time, 1-frequency_dot_2(peak_x), 'x'); % 标注峰值
            ylabel('波形变化');
            title('请调整海明窗宽度 STFT window width, 重叠长度 STFT noverlap, 以提高时域分辨率');
            
            subplot(4,1,3);
            plot((1:length(obj.clips))/obj.sample_rate, obj.clips);
            hold on;
            for i = obj.divide_time
                plot([i, i], [-0.5, 0.5], color = 'r');
            end
            obj.divide_time = [obj.divide_time; length(obj.clips)/obj.sample_rate];
            obj.middle_time = (1-0.4)*obj.divide_time(1:end-1) + 0.4*obj.divide_time(2:end);
            for i = obj.middle_time
                plot([i - obj.T_analyse_length/2/obj.sample_rate, ...
                      i + obj.T_analyse_length/2/obj.sample_rate, ...
                      i + obj.T_analyse_length/2/obj.sample_rate, ...
                      i - obj.T_analyse_length/2/obj.sample_rate, ...
                      i - obj.T_analyse_length/2/obj.sample_rate, ...
                      ], [-0.1, -0.1, 0.1, 0.1, -0.1], color = 'r');
            end
            xlim([0, length(obj.clips)/obj.sample_rate]);
            title('请调整取内积间隔 dot step, 最小峰值 dot min peak height, 最小峰间距 dot min peak distance, 以避免漏标记');
        end

        function [frequency, harmonic] = frequency_harmonic(obj, single_tone_clips)
            % 音调、谐波分析

            % 估计基音周期
            R_single_tone = xcorr(single_tone_clips); % 自相关函数
            R_single_tone = xcorr(R_single_tone);
            R_single_tone_length = length(R_single_tone);

            x_shift = obj.sample_rate / R_single_tone_length * ...
                (-floor(R_single_tone_length/2) : ceil(R_single_tone_length/2-1))';
            DFT = abs(fftshift(fft(R_single_tone)));
            try
                [peak_y, peak_x] = findpeaks(DFT);
                base_point = peak_x(peak_y == max(peak_y));
                approx_frequency = x_shift(base_point(2));
                approx_T = obj.sample_rate / approx_frequency;
            catch
                error('估计基音周期失败');
            end
            if approx_frequency < music_const.frequency.f
                approx_T = approx_T / 2;
                warning('估计基音周期可能不准确, 音调估计值过低');
            elseif approx_frequency > music_const.frequency.e*3
                approx_T = approx_T * 2;
                warning('估计基音周期可能不准确, 音调估计值过高');
            end

            % 确定基音周期
            resample_wave = resample(single_tone_clips, obj.resample_times, 1);
            R_resample_wave = xcorr(resample_wave);
            [~, peak_x] = findpeaks(R_resample_wave, ...
                MinPeakDistance = approx_T * obj.resample_times * obj.T_time_threshold);
        
            point_count   = floor(length(peak_x)/2);
            point_distance = ceil(length(peak_x)/2);
            resample_wave_T = (sum(peak_x(end-point_count+1:end)) - sum(peak_x(1:point_count))) ...
                / point_distance / point_count; % 重采样波形一周期样点数
        
            % 时间同步平均
            resample_wave_length = length(resample_wave); % 重采样波形样点数
            resample_wave_single_T = tsa(resample_wave, 1, 0:resample_wave_T:resample_wave_length);

            % 频率
            frequency = 1 / (resample_wave_T / obj.resample_times / obj.sample_rate);

            % 周期重复 DFT
            repeat_resample_wave = repmat(resample_wave_single_T, obj.DFT_repeat_times, 1);
            repeat_resample_wave = repeat_resample_wave / max(resample_wave_length);
            repeat_resample_wave_length = length(repeat_resample_wave);
            x_shift = obj.sample_rate * obj.DFT_repeat_times / repeat_resample_wave_length * ...
                (-floor(repeat_resample_wave_length/2) : ceil(repeat_resample_wave_length/2-1))';
            DFT = abs(fftshift(fft(repeat_resample_wave)));
            base_frequency_peak = max(DFT(x_shift>frequency-obj.DFT_base_frequency_range & ...
                x_shift<frequency+obj.DFT_base_frequency_range));
            if base_frequency_peak < max(DFT)/10
                warning('基频幅值过低, 确定基频可能出错, 请增大 DFT_base_frequency_range 并重试');
            end

            % 谐波分量
            harmonic = findpeaks(DFT(x_shift>frequency+obj.DFT_base_frequency_range & ...
                x_shift<obj.DFT_max_frequency) / base_frequency_peak, ...
                MinPeakHeight=obj.DFT_min_harmonic);
        end

        function beat_analyse(obj)
            % 确定节拍

            % 合并重复的分隔
            last_tone = obj.tone{1};
            for i = 2 : length(obj.tone)
                if strcmp(obj.tone{i}, last_tone)
                    if max(abs(obj.clips(round(obj.divide_time(i)*obj.sample_rate-obj.T_analyse_length) : ...
                                       round(obj.divide_time(i)*obj.sample_rate)))) ...
                     > obj.beat_volume_threshold_next * max(abs(obj.clips(round(obj.divide_time(i)*obj.sample_rate) : ...
                                       round(obj.divide_time(i)*obj.sample_rate+obj.T_analyse_length))))
                        % 基频不变, 并且音量不突然增大
                        obj.divide_time(i) = 0;
                    end
                else
                    last_tone = obj.tone{i};
                end
            end
            harmonic_struct = struct();
            for i = 1:length(obj.tone)
                if isfield(harmonic_struct, obj.tone{i})
                    length_main = length(harmonic_struct.(obj.tone{i}));
                    length_add = length(obj.harmonic{i});
                    difference = length_main - length_add;
                    if difference > 0
                        harmonic_struct.(obj.tone{i}) = harmonic_struct.(obj.tone{i}) + ...
                            [obj.harmonic{i}; zeros(difference, 1)];
                    elseif difference < 0
                        harmonic_struct.(obj.tone{i}) = obj.harmonic{i} + ...
                            [harmonic_struct.(obj.tone{i}); zeros(-difference, 1)];
                    else
                        harmonic_struct.(obj.tone{i}) = harmonic_struct.(obj.tone{i}) + obj.harmonic{i};
                    end
                else
                    harmonic_struct.(obj.tone{i}) = obj.harmonic{i};
                end
            end
            obj.tone = obj.tone(obj.divide_time(1:end-1) > 0);
            for i = 1:length(obj.tone)
                harmonic_array = harmonic_struct.(obj.tone{i});
                harmonic_array = harmonic_array / harmonic_array(1);
                harmonic_struct.(obj.tone{i}) = harmonic_array(2:end);
            end
            obj.harmonic = harmonic_struct;
            obj.divide_time = obj.divide_time(obj.divide_time > 0);
            subplot(4,1,4);
            plot((1:length(obj.clips))/obj.sample_rate, obj.clips);
            hold on;
            for i = obj.divide_time
                plot([i, i], [-0.5, 0.5], color = 'r');
            end
            xlabel('时间 / s');
            xlim([0, length(obj.clips)/obj.sample_rate]);
            title('合并部分分隔点');
            
            divide_length = obj.divide_time(2:end) - obj.divide_time(1:end-1);
            reference_divide_length = median(divide_length);
            obj.beat = (divide_length / reference_divide_length);
            while true
                if max(round(obj.beat)) > 4
                    obj.beat = (obj.beat / 2);
                else
                    obj.beat(obj.beat>=1) = round(obj.beat(obj.beat>=1));
                    obj.beat(obj.beat<1) = 2.^round(log2(obj.beat(obj.beat<1)));
                    break;
                end
            end
        end
    end

    methods (Static)
        function tone = tone_recognition(frequency)
            if frequency*2 < music_const.frequency.B + music_const.frequency.c
                % 大字组
                distance = floor(12 * log2(frequency/music_const.frequency.B));
                if frequency*2 > music_const.frequency.B*(2^(distance/12)+2^((distance+1)/12))
                    distance = distance + 1; % 频率线性距离接近后一音
                end
                group = floor(-distance / 12);
                letter = music_const.Tone_Letter{mod(-distance, 12)+1};
            else
                % 小字组
                distance = floor(12 * log2(frequency/music_const.frequency.c));
                if frequency*2 > music_const.frequency.c*(2^(distance/12)+2^((distance+1)/12))
                    distance = distance + 1; % 频率线性距离接近后一音
                end
                group = floor(distance / 12);
                letter = music_const.tone_letter{mod(distance, 12)+1};
            end
            tone = [letter, num2str(group)];
        end
    end

    methods (Static, Access = private)
        function set_callback(~, event)
            try
                event.AffectedObject.analyse();
            catch
                warning('music_analysis 更新失败');
            end
        end
    end
    
end