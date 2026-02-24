classdef FxLMS < AdaptiveFilter
    properties
        eps = 1e-8;  % small constant to avoid division by zero
        step_size = 0.05;
    end

    methods
        function obj = FxLMS(ref, target, ir, filter_order, step_size, varargin)
            obj@AdaptiveFilter(ref, target, ir, filter_order, 'FxLMS');
            
            obj.step_size = step_size;
            if nargin > 6
                obj.eps = varargin{1};
            end
        end

        
        function to_gpu(obj); end
        function to_cpu(obj); end

        
        function update_coefficients(obj)
            x = zeros(obj.filter_order, 1);

            if obj.gpu
                x = gpuArray(x);               
            end

            xnorm = 0;
        
            for n = 1:length(obj.target)

                old = x(end);

                x(2:end) = x(1:end-1);
                x(1) = obj.x_filt(n);

                xnorm = xnorm + x(1)^2 - old^2;

                obj.y(n) = obj.w.' * x;
                obj.e(n) = obj.target(n) - obj.y(n);
            
                obj.w = obj.w + (obj.step_size / (obj.eps + xnorm)) * obj.e(n) * x;
            end
        end

    end
end