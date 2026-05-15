classdef FxLMS < AdaptiveFilter
    properties
        eps = 1e-8;  % small constant to avoid division by zero
        step_fn;
        use_norm = true;

        x;
        xnorm = 0;
    end

    methods
        function obj = FxLMS(ref, target, ir, filter_order, step_size, use_norm, varargin)
            obj@AdaptiveFilter(ref, target, ir, filter_order, 'FxLMS');
            
            obj.step_fn = step_size;
            obj.use_norm = use_norm;

            if nargin > 6
                obj.eps = varargin{1};
            end

            obj.x = zeros(obj.filter_order, 1);
            obj.xnorm = 0;
        end

        function val = step_size(obj, varargin)
            if isnumeric(obj.step_fn)
                val = obj.step_fn;
                return;
            end
            val = obj.step_fn(varargin{1});

        end

        
        function to_gpu(obj); end
        function to_cpu(obj); end

        
        function update_coefficients(obj, iter)
            step = obj.step_size(iter);

            for n = 1:length(obj.target)

                old = obj.x(end);

                obj.x(2:end) = obj.x(1:end-1);
                obj.x(1) = obj.x_filt(n);


                obj.y(n) = obj.w.' * obj.x;
                obj.e(n) = obj.target(n) - obj.y(n);
                
                if obj.use_norm
                    obj.xnorm = obj.xnorm + obj.x(1)^2 - old^2;
                    obj.w = obj.w + (step / (obj.eps + obj.xnorm)) * obj.e(n) * obj.x;
                else
                    obj.w = obj.w + step * obj.e(n) * obj.x;
                end

            end
        end

    end
end