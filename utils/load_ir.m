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