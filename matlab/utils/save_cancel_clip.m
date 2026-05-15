function save_cancel_clip(wav, filename)
    [a, fs] = arthur();
    
    wav = wav(1:length(a));
    y = ([wav, a]);
    no_panel = [zeros(length(wav), 1), a];

    cancel = ['audio_experiments', filesep, filename, '_cancel.wav'];
    nc = ['audio_experiments', filesep, filename, '_nc.wav'];
    audiowrite(cancel, y, fs);
    audiowrite(nc, no_panel, fs);
end