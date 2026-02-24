classdef FxRLS < AdaptiveFilter
    properties
        lamda = 0.999; % larger = longer memory over signal
        cov_scale = 0.1; % larger = more aggressive convergence
        P;

        y;
        e;
    end

    methods
        function obj = FxRLS(ref, target, ir, filter_order, lamda, cov_scale)
            obj@AdaptiveFilter(ref, target, ir, filter_order, 'FxRLS');

            obj.lamda = lamda;
            obj.cov_scale = cov_scale;
            obj.P = eye(filter_order) * cov_scale;

            
            
        end

        function reset(obj)
            reset@AdaptiveFilter(obj);
            obj.P = eye(obj.filter_order) * obj.cov_scale;
        end

        function to_gpu(obj)
            obj.P = gpuArray(obj.P);
        end
        
        function to_cpu(obj)
            obj.P = gather(obj.P);
        end

        
        function update_coefficients(obj)
            x = zeros(obj.filter_order, 1);

            if obj.gpu
                x = gpuArray(x);            
            end
            
            for n = 1:length(obj.target)
                x(2:end) = x(1:end-1);
                x(1) = obj.x_filt(n);

                obj.y(n) = obj.w.' * x;
                obj.e(n) = obj.target(n) - obj.y(n);
            
                % RLS update
                Px = obj.P * x;
                denom = obj.lamda + x.' * Px;     
                k = Px / denom; % gain vector (N x 1)
    
                obj.w = obj.w + obj.k * obj.e(n);
                obj.P = (obj.P - k * (x.' * obj.P)) / obj.lamda; % covariance update
            end
        end
    end
end