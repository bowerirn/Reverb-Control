function [reverb, ref, fs] = recorded_samples(~)
    if nargin == 0
        [mics, fs] = audioread("audio_experiments/mics_nc.wav");
    else
        [mics, fs] = audioread("audio_experiments/mics_foam.wav");
    end
    ref = mics(:, 2);
    reverb = mics(:, 1);
end