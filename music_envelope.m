function envelope = music_envelope(time_split, peak, change_rate, sample_rate)
% 参数 time_split = [t1, t2, t3, t4, t5] 为 5 个阶段的分界点
% 参数 peak 为冲激峰值
% 参数 change_rate = [k1, k2, k3, k4] 控制指数变化速率，其中 k1, k3 > 0; k2, k4 < 0
% 参数 sample_rate 为采样率

time_split = round(time_split * sample_rate);
change_rate = change_rate / sample_rate;

time_period_1 = (0:time_split(1)-1)';
time_period_2 = (time_split(1):time_split(2)-1)';
time_period_4 = (time_split(3):time_split(4)-1)';
time_period_5 = (time_split(4):time_split(5)-1)';

envelope = [
    peak * (exp(change_rate(1)*time_period_1)-1) ./ (exp(change_rate(1)*time_split(1))-1);
    
    peak * (exp(change_rate(2)*time_period_2)-exp(change_rate(2)*time_split(2))) ...
        ./ (exp(change_rate(2)*time_split(1))-exp(change_rate(2)*time_split(2))) ...
    + 1 *  (exp(change_rate(2)*time_split(1))-exp(change_rate(2)*time_period_2)) ...
        ./ (exp(change_rate(2)*time_split(1))-exp(change_rate(2)*time_split(2)));
    
    ones(time_split(3)-time_split(2),1);

    log(2) * (exp(change_rate(3)*time_split(3))-exp(change_rate(3)*time_period_4)) ...
           ./ (exp(change_rate(3)*time_split(3))-exp(change_rate(3)*time_split(4))) ...
    + 1 *     (exp(change_rate(3)*time_period_4)-exp(change_rate(3)*time_split(4))) ...
           ./ (exp(change_rate(3)*time_split(3))-exp(change_rate(3)*time_split(4)));
    
    log(2) * (exp(change_rate(4)*time_period_5)-exp(change_rate(4)*time_split(5))) ...
           ./ (exp(change_rate(4)*time_split(4))-exp(change_rate(4)*time_split(5)));
    ] / peak;

end

