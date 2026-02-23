function [wav, fs] = arthur()
    wav, fs = audioread('data/arthur_clip_48k.wav');
end

function [reverb, ref, fs] = recorded_data()
    [mics, fs] = audioread("double_mics.wav");
    ref = mics(:, 2);
    reverb = mics(:, 1);
end


function [panel, ref, reverb, time] = irs()
    load('data/ir_delay_tests.mat')
    panel = measurementData(1, 2).ImpulseResponse.Amplitude;
    ref = measurementData(2, 2).ImpulseResponse.Amplitude;
    reverb = measurementData(3, 2).ImpulseResponse.Amplitude;
    time = measurementData(1, 2).ImpulseResponse.Time;
end

function [ir, time] = load_ir(ir_name)
        load('data/ir_delay_tests.mat')
        switch ir_name
            case "panel"
                ir = measurementData(1, 2).ImpulseResponse.Amplitude;
                time = measurementData(1, 2).ImpulseResponse.Time;
            case "ref"
                ir = measurementData(2, 2).ImpulseResponse.Amplitude;
                time = measurementData(2, 2).ImpulseResponse.Time;
            case "reverb"
                ir = measurementData(3, 2).ImpulseResponse.Amplitude;
                time = measurementData(3, 2).ImpulseResponse.Time;
            otherwise
                error("Unknown IR name: %s", ir_name);
        end
end