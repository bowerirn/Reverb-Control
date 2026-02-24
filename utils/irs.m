function [panel, ref, reverb, fs, time] = irs()
    load('data/ir_delay_tests.mat');
    panel = measurementData(1, 2).ImpulseResponse.Amplitude;
    ref = measurementData(2, 2).ImpulseResponse.Amplitude;
    reverb = measurementData(3, 2).ImpulseResponse.Amplitude;

    fs = measurementData(1, 4).SampleRate;
    time = measurementData(1, 2).ImpulseResponse.Time;
end