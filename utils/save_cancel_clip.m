function save_cancel_clip(wav, filename)
    [a, fs] = arthur();
    
        disp(size(a));
        disp(size(wav));

    wav = wav(1:length(a));

    y = ([wav, a]);
    audiowrite([filename, '.wav'], y, fs);
    audiowrite([filename, '_no_cancel.wav'], [zeros(size(wav), 1), a], fs);
end