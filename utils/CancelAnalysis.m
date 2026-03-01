classdef CancelAnalysis < handle
    properties
        name;
        error_mic;
        error_mic_nc
        source_mic = false;
        source_mic_nc = false;
        fs;

        cancel_rms;
        nc_rms;

        t;

    end


    methods
        function obj = CancelAnalysis(name)
            obj.name = name;
            cancel_fn = ['audio_experiments', filesep, name, '_mics.wav'];
            nc_fn = ['audio_experiments', filesep, name, '_mics_nc.wav'];
        
            [mics, fs] = audioread(cancel_fn);
            obj.error_mic = mics(:, 1);
            obj.fs = fs;
        
            if size(mics, 2) == 2
                obj.source_mic = mics(:, 2);
            end


            [mics_nc, fs_nc] = audioread(nc_fn);
            obj.error_mic_nc = mics_nc(:, 1);
        
            if size(mics_nc, 2) == 2
                obj.source_mic_nc = mics_nc(:, 2);
            end

            obj.cancel_rms = rms(obj.error_mic);
            obj.nc_rms = rms(obj.error_mic_nc);

            obj.t = (0:length(obj.error_mic) - 1) / fs;
        end

        function plot_cancel_results(obj)
            figure;

            if obj.cancel_rms > obj.nc_rms
                plot(obj.t, obj.error_mic);
                hold on;
                plot(obj.t, obj.error_mic_nc);
                hold off;
                legend('With Cancel', 'No Cancel');
            else
                plot(obj.t, obj.error_mic_nc);
                hold on;
                plot(obj.t, obj.error_mic);
                hold off;
                legend('No Cancel', 'With Cancel');
            end

            title(['Cancellation Results for ', obj.name]);
            xlabel('Time (s)');
            ylabel('Amplitude');
            xlim([obj.t(1), obj.t(end)]);
        end
    end




end