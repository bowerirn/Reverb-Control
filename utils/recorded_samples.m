function [reverb, ref, fs] = recorded_samples()
    [mics, fs] = audioread("double_mics.wav");
    ref = mics(:, 2);
    reverb = mics(:, 1);
end