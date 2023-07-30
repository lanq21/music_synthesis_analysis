classdef music_synthesis < handle
    
    properties (SetObservable, AbortSet)
        sample_rate = 8000                                % 采样率
        speed = music_const.speed.('Moderato')            % 速度(BMP)
        beat_tempo                                        % 节拍
        base_frequency                                    % 唱名 "1" 频率
        envelope_description = music_const.envelope_piano % 包络描述
        harmonic = music_const.harmonic_piano             % 谐波分量
        clips_description = '1 2 3'                       % 片段简谱
    end

    properties (Dependent)
        full_duration    % 音符时长, 不包含迭接部分
        tail_length      % 迭接部分长度
        beat             % 每小节有 beat 拍
        tempo            % tempo 分音符为一拍
    end

    properties (SetAccess = private)
        time             % {全音符, 1/2 音符, 1/4 音符, 1/8 音符, 1/16 音符, 1/32 音符} 时间序列, 包含迭接部分
        time_dot         % 附点 {全音符, 1/2 音符, 1/4 音符, 1/8 音符, 1/16 音符, 1/32 音符} 时间序列, 包含迭接部分
        frequency        % (唱名).normal / flat / sharp 频率
        envelope         % {全音符, 1/2 音符, 1/4 音符, 1/8 音符, 1/16 音符, 1/32 音符} 包络
        envelope_dot     % 附点 {全音符, 1/2 音符, 1/4 音符, 1/8 音符, 1/16 音符, 1/32 音符} 包络
        tone             % .t_分音符.f_音调 单音波形
        tone_dot         % .t_分音符.f_音调 附点单音波形
        clips = []       % 片段波形
    end

    properties (Access = private)
        time_up_to_date = 0
        time_dot_up_to_date = 0
        frequency_up_to_date = 0
        envelope_up_to_date = 0
        envelope_dot_up_to_date = 0
        tone_up_to_date = 0
        tone_dot_up_to_date = 0
        clips_up_to_date = 0
    end

    methods % 初始化方法
        function obj = music_synthesis(varargin)
            n = nargin;
            for i = 1:2:n
                varable_name_cell  = varargin(i);
                varable_name       = varable_name_cell{1};
                varable_cell       = varargin(i+1);
                varable            = varable_cell{1};
                obj.(varable_name) = varable;
            end
            addlistener(obj, 'sample_rate',          'PostSet', @music_synthesis.set_callback);
            addlistener(obj, 'speed',                'PostSet', @music_synthesis.set_callback);
            addlistener(obj, 'beat_tempo',           'PostSet', @music_synthesis.set_callback);
            addlistener(obj, 'base_frequency',       'PostSet', @music_synthesis.set_callback);
            addlistener(obj, 'envelope_description', 'PostSet', @music_synthesis.set_callback);
            addlistener(obj, 'harmonic',             'PostSet', @music_synthesis.set_callback);
            addlistener(obj, 'clips_description',    'PostSet', @music_synthesis.set_callback);
        end
    end

    methods % SetObservable 属性 set 方法
        function set.speed(obj, speed_description)
            if isscalar(speed_description)
                obj.speed = speed_description;
            else
                obj.speed = music_const.speed.(speed_description);
            end
        end

        function set.base_frequency(obj, tune_description)
            modulation = 1; % 变调标记: 无升调 / 降调
            group = 0;      % 组: 小字组 / 大字组
            if strlength(tune_description) == 3      % 有变调、组记号
                if strcmp(tune_description(1), '#')
                    modulation = 2^(1/12);
                else
                    modulation = 2^(-1/12);
                end
                freq = music_const.frequency.(tune_description(2));
                group = str2double(tune_description(3));
            elseif strlength (tune_description) == 1 % 无记号
                freq = music_const.frequency.(tune_description);
            else
                if isletter(tune_description(2)) % 有变调记号
                    if strcmp(tune_description(1), '#')
                        modulation = 2^(1/12);
                    else
                        modulation = 2^(-1/12);
                    end
                    freq = music_const.frequency.(tune_description(2));
                else                           % 有组记号
                    freq = music_const.frequency.(tune_description(1));
                    group = str2double(tune_description(2));
                end
            end
            if freq < 127 % 大字 group 组
                obj.base_frequency = modulation * freq * 2^(-group);
            else          % 小字 group 组
                obj.base_frequency = modulation * freq * 2^(group);
            end
        end
    end

    methods % Dependent 属性 get 方法
        function value = get.full_duration(obj)
            value = 60 / obj.speed * obj.tempo;
        end

        function value = get.beat(obj)
            value = str2double(obj.beat_tempo(1));
        end

        function value = get.tempo(obj)
            value = str2double(obj.beat_tempo(3));
        end
    end

    methods % private 属性 get 方法
        function value = get.time(obj)
            if obj.time_up_to_date == 0
                for i = 2.^(0:5)
                    len = length(obj.envelope{i});
                    obj.time{i} = (0 : 1/obj.sample_rate : (len-1) / obj.sample_rate)';
                end
            end
            value = obj.time;
        end

        function value = get.time_dot(obj)
            if obj.time_dot_up_to_date == 0
                for i = 2.^(0:5)
                    len = length(obj.envelope_dot{i});
                    obj.time_dot{i} = (0 : 1/obj.sample_rate : (len-1) / obj.sample_rate)';
                end
            end
            value = obj.time_dot;
        end

        function value = get.frequency(obj)
            if obj.frequency_up_to_date == 0
                obj.frequency = struct( ...
                    'normal', num2cell(obj.base_frequency * 2.^([0,2,4,5,7,9,11]/12)), ...
                    'flat',   num2cell(obj.base_frequency * 2.^([0,2,4,5,7,9,11]/12 - 1)), ...
                    'sharp',  num2cell(obj.base_frequency * 2.^([0,2,4,5,7,9,11]/12 + 1))  ...
                    );
                obj.frequency_up_to_date = 1;
            end
            value = obj.frequency;            
        end

        function value = get.envelope(obj)
            if obj.envelope_up_to_date == 0
                for i = 2.^(0:5)
                    obj.envelope{i} = music_envelope( ...
                        obj.envelope_description(1:5) * obj.full_duration/i, ...
                        obj.envelope_description(6), ...
                        obj.envelope_description(7:10) * i, ...
                        obj.sample_rate ...
                    );
                    obj.envelope{64} = music_envelope( ...
                        obj.envelope_description(1:5) * obj.full_duration/obj.tempo * 1.5, ...
                        obj.envelope_description(6), ...
                        obj.envelope_description(7:10) * obj.tempo / 1.5, ...
                        obj.sample_rate ...
                    );
                end
                obj.envelope_up_to_date = 1;
            end
            value = obj.envelope;
        end

        function value = get.envelope_dot(obj)
            if obj.envelope_dot_up_to_date == 0
                for i = 2.^(0:5)
                    obj.envelope_dot{i} = music_envelope( ...
                        obj.envelope_description(1:5) * obj.full_duration/i * 1.5, ...
                        obj.envelope_description(6), ...
                        obj.envelope_description(7:10) * i / 1.5, ...
                        obj.sample_rate ...
                    );
                end
                obj.envelope_dot_up_to_date = 1;
            end
            value = obj.envelope_dot;
        end

        function value = get.tone(obj)
            if obj.tone_up_to_date == 0
                for i = 2.^(0:5) % 时长
                    for j = 1:7  % 音调
                        tone_normal = sin(2*pi*obj.frequency(j).normal * obj.time{i});
                        tone_flat   = sin(2*pi*obj.frequency(j).flat   * obj.time{i});
                        tone_sharp  = sin(2*pi*obj.frequency(j).sharp  * obj.time{i});
                        for k = 1 : length(obj.harmonic) % 谐波
                            tone_normal = tone_normal + ...
                                obj.harmonic(k) * sin(2*pi*obj.frequency(j).normal * (k+1) * obj.time{i});
                            tone_flat   = tone_flat   + ...
                                obj.harmonic(k) * sin(2*pi*obj.frequency(j).flat   * (k+1) * obj.time{i});
                            tone_sharp  = tone_sharp  + ...
                                obj.harmonic(k) * sin(2*pi*obj.frequency(j).sharp  * (k+1) * obj.time{i});
                        end
                        obj.tone.(['t_', num2str(i)]).('f_0') = ...
                            zeros(size(obj.time{i}));
                        obj.tone.(['t_', num2str(i)]).(['f_', num2str(j)     ]) = ...
                            obj.envelope{i} .* tone_normal / max(abs(tone_normal));
                        obj.tone.(['t_', num2str(i)]).(['f_', num2str(j), 'f']) = ...
                            obj.envelope{i} .* tone_flat   / max(abs(tone_flat  ));
                        obj.tone.(['t_', num2str(i)]).(['f_', num2str(j), 's']) = ...
                            obj.envelope{i} .* tone_sharp  / max(abs(tone_sharp ));
                    end
                end
                obj.tone_up_to_date = 1;
            end
            value = obj.tone;
        end

        function value = get.tone_dot(obj)
            if obj.tone_dot_up_to_date == 0
                for i = 2.^(0:5) % 时长
                    for j = 1:7  % 音调
                        tone_normal = sin(2*pi*obj.frequency(j).normal * obj.time_dot{i});
                        tone_flat   = sin(2*pi*obj.frequency(j).flat   * obj.time_dot{i});
                        tone_sharp  = sin(2*pi*obj.frequency(j).sharp  * obj.time_dot{i});
                        for k = 1 : length(obj.harmonic) % 谐波
                            tone_normal = tone_normal + ...
                                obj.harmonic(k) * sin(2*pi*obj.frequency(j).normal * (k+1) * obj.time_dot{i});
                            tone_flat   = tone_flat   + ...
                                obj.harmonic(k) * sin(2*pi*obj.frequency(j).flat   * (k+1) * obj.time_dot{i});
                            tone_sharp  = tone_sharp  + ...
                                obj.harmonic(k) * sin(2*pi*obj.frequency(j).sharp  * (k+1) * obj.time_dot{i});
                        end
                        obj.tone_dot.(['t_', num2str(i)]).('f_0') = ...
                            zeros(size(obj.time{i}));
                        obj.tone_dot.(['t_', num2str(i)]).(['f_', num2str(j)     ]) = ...
                            obj.envelope_dot{i} .* tone_normal / max(abs(tone_normal));
                        obj.tone_dot.(['t_', num2str(i)]).(['f_', num2str(j), 'f']) = ...
                            obj.envelope_dot{i} .* tone_flat   / max(abs(tone_flat  ));
                        obj.tone_dot.(['t_', num2str(i)]).(['f_', num2str(j), 's']) = ...
                            obj.envelope_dot{i} .* tone_sharp  / max(abs(tone_sharp ));
                    end
                end
                obj.tone_dot_up_to_date = 1;
            end
            value = obj.tone_dot;
        end

        function value = get.clips(obj)
            if obj.clips_up_to_date == 0
                clips_cell = strsplit(strip(obj.clips_description), ' ');
                tone_count   = 1; % '1-' 作为一个音节计数
                len = length(clips_cell);
                tone_t = zeros(len, 1);
                tone_f = cell(len, 1);
                dot = zeros(len, 1); % 附点音符标记
                delay_to_dot = 0; % 延时到附点音符
                for i = 1 : len
                    if strcmp(clips_cell{i}, '-')
                        if delay_to_dot
                            dot(tone_count - 1) = 1;
                        else
                            dot(tone_count - 1) = 0;
                            tone_t(tone_count-1) = tone_t(tone_count-1) / 2;
                        end
                        delay_to_dot = ~delay_to_dot;
                    elseif strcmp(clips_cell{i}, '.')
                        delay_to_dot = 0;
                        dot(tone_count - 1) = 1;
                    else
                        delay_to_dot = 0;
                        tone_t(tone_count) = obj.tempo;
                        for j = 1:5
                            if strcmp(clips_cell{i}(j), '_')
                                tone_t(tone_count) = tone_t(tone_count) * 2;
                            else
                                break;
                            end
                        end
                        if strcmp(clips_cell{i}(j), 'b') || strcmp(clips_cell{i}(j), 'f')
                            tone_f{tone_count} = ['f_', clips_cell{i}(j+1), 'f'];
                        elseif strcmp(clips_cell{i}(j), '#') || strcmp(clips_cell{i}(j), 's')
                            tone_f{tone_count} = ['f_', clips_cell{i}(j+1), 's'];
                        else
                            tone_f{tone_count} = ['f_', clips_cell{i}(j)];
                        end
                        tone_count = tone_count + 1;
                    end
                end

                obj.clips = [];
                tone_tail = [];
                for i = 1 : tone_count-1
                    if dot(i)
                        tone_i = obj.tone_dot.(['t_', num2str(tone_t(i))]).(tone_f{i});
                        main_tail_cut_off = round(obj.full_duration * obj.sample_rate / tone_t(i) * 1.5);
                    else
                        tone_i = obj.tone.(['t_', num2str(tone_t(i))]).(tone_f{i});
                        main_tail_cut_off = round(obj.full_duration * obj.sample_rate / tone_t(i));
                    end
                    
                    tone_main = tone_i(1 : main_tail_cut_off - 1);
                    obj.clips = [ ...
                        obj.clips; ...
                        [tone_tail; zeros(length(tone_main)-length(tone_tail), 1)] + tone_main];
                    tone_tail = tone_i(main_tail_cut_off : end);
                end
                obj.clips = [obj.clips; tone_tail];

                obj.clips_up_to_date = 1;
            end
            obj.clips = obj.clips / max(abs(obj.clips));
            value = obj.clips;
        end
    end

    methods (Static, Access = private)
        function set_callback(src, event)
            event.AffectedObject.clips_up_to_date = 0;

            if strcmp(src.Name, 'clips_description') == 0
                event.AffectedObject.tone_up_to_date = 0;
                event.AffectedObject.tone_dot_up_to_date = 0;
            end
            if strcmp(src.Name, 'envelope_description')
                event.AffectedObject.envelope_up_to_date = 0;
                event.AffectedObject.envelope_dot_up_to_date = 0;
                event.AffectedObject.time_up_to_date = 0;
                event.AffectedObject.time_dot_up_to_date = 0;
            elseif strcmp(src.Name, 'base_frequency')
                event.AffectedObject.frequency_up_to_date = 0;
            elseif strcmp(src.Name, 'sample_rate')
                event.AffectedObject.time_up_to_date = 0;
                event.AffectedObject.time_dot_up_to_date = 0;
            end
        end
    end

end