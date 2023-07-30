function clips = generate_clips(tones, table, sample_rate, overlapping_length)
    tail_duration = overlapping_length * sample_rate;
    clips = [];
    tail = zeros(tail_duration,1);
    for i = 1:length(table)   
        tone = tones{28 + table(i)};
        clips = [ ...
            clips; ...
            [tail;zeros(length(tone)-2*tail_duration,1)]+tone(1:end-tail_duration) ...
            ];
        tail = tone(end-tail_duration+1:end);    
    end
    clips = [clips;tail];
end