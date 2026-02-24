classdef (Abstract) AdaptiveFilter < handle
    properties
        ref
        target
        ir
        filter_order
        name

        gpu = gpuDeviceCount > 0;

        y;
        e;

        w
        x_filt
        ref_filt
    end
    
    methods
        function obj = AdaptiveFilter(ref, target, ir, filter_order, name)
            obj.ref = ref;
            obj.target = target;
            obj.ir = ir;
            obj.filter_order = filter_order;
            obj.name = name;

            obj.w = zeros(filter_order, 1);
            obj.x_filt = filter(ir, 1, ref);

            obj.y = zeros(1, length(obj.target));
            obj.e = zeros(1, length(obj.target));
        end

        function value = get.ref_filt(obj)
             value = filter(obj.w, 1, obj.ref);
        end

        function reset(obj)
            obj.w = zeros(obj.filter_order, 1);
        end

        function results_plot(obj)
            pred = filter(obj.w, 1, obj.x_filt);
            figure
            plot(obj.target)
            hold on
            plot(pred)
            legend('ground truth', 'learned prediction')
            title(obj.name)
        end
        

        function [final_pred, mses] = learn(obj, n_iter, varargin)
            log = nargin > 2;

            mses = zeros(1, n_iter);

            if obj.gpu
                obj.w = gpuArray(obj.w);
                obj.target = gpuArray(obj.target);
                obj.ref = gpuArray(obj.ref);
                obj.x_filt = gpuArray(obj.x_filt);
                obj.y = gpuArray(obj.y);
                obj.e = gpuArray(obj.e);
                obj.to_gpu();
            end

            for i = 1:n_iter
                if log
                    fprintf("Iteration %d/%d\n", i, n_iter);
                end
                obj.update_coefficients();
                mses(i) = mean(gather(obj.e).^2);
            end

            if obj.gpu
                obj.w = gather(obj.w);
                obj.target = gather(obj.target);
                obj.ref = gather(obj.ref);
                obj.x_filt = gather(obj.x_filt);
                obj.y = gather(obj.y);
                obj.e = gather(obj.e);
                obj.to_cpu();
            end

            final_pred = obj.y;
        end
        
        
    end

    methods (Abstract)
        [y, e] = update_coefficients(obj);
        to_gpu(obj);
        to_cpu(obj);

    end
end