classdef (Abstract) AdaptiveFilter < handle
    properties
        ref
        target
        ir
        filter_order
        step_size

        w
        x_filt
        ref_filt
    end
    
    methods
        function obj = AdaptiveFilter(ref, target, ir, filter_order, step_size)
            obj.ref = ref;
            obj.target = target;
            obj.ir = ir;
            obj.filter_order = filter_order;
            obj.step_size = step_size;

            obj.w = zeros(filter_order, 1);
            obj.x_filt = filter(ir, 1, ref);
        end

        function value = get.ref_filt(obj)
             value = filter(obj.w, 1, obj.ref);
        end
        
        function output = filter(obj, signal)
            output = filter(obj.w, 1, signal);
        end

        function [final_pred, mses] = learn(obj, n_iter, varargin)
            log = nargin > 2;

            mses = zeros(1, n_iter);

            for i = 1:n_iter
                if log
                    fprintf("Iteration %d/%d\n", i, n_iter);
                end
                [y, e] = obj.update_coefficients();
                mses(i) = mean(e.^2);
            end
        end
        
        
    end

    methods (Abstract)
        [y, e] = update_coefficients(obj);
    end
end