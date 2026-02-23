classdef FxLMS < AdaptiveFilter
    properties
        eps = 1e-8;  % small constant to avoid division by zero
    end

    methods
        function obj = FxLMS(ref, target, ir, filter_order, step_size, varargin)
            obj@AdaptiveFilter(ref, target, ir, filter_order, step_size);

            if nargin > 6
                obj.eps = varargin{1};
            end
        end
        
        function [y, e] = update_coefficients(obj)
            x = zeros(filter_order, 1);
            y = zeros(1, length(obj.target));
            e = zeros(1, length(target));
            xnorm = 0;
        
            for n = 1:length(target)

                old = x(end);

                x(2:end) = x(1:end-1);
                x(1) = obj.x_filt(n);

                xnorm = xnorm + x(1)^2 - old^2;

                y(n) = obj.w.' * obj.x;
                e(n) = obj.target(n) - y(n);
            
                obj.w = obj.w + (obj.step_size / (obj.eps + xnorm)) * e(n) *x;
            end
        end
    end
end