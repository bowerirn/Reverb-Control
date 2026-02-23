function [u, w, y, mse] = fxlms(ref, target, ir, n_iter, filter_order, mu, eps)

    w = zeros(filter_order, 1);
    y = zeros(1, length(target)); % the final output of (x_ref * w)[n]
    e = zeros(1, length(target)); % the error 
    
    x_filt = filter(ir, 1, ref);
    
    mse = zeros(n_iter, 1);
    
    for i = 1:n_iter
    
        x = zeros(filter_order, 1);
        xnorm = 0;
    
        for n = 1:length(target)

            old = x(end);

            x(2:end) = x(1:end-1);
            x(1) = x_filt(n);

            xnorm = xnorm + x(1)^2 - old^2; % faster than computing the full dot product every time

            y(n) = w.' * x;
            e(n) = target(n) - y(n);
           
            w = w + (mu / (eps + xnorm)) * e(n) *x;
        end
    
        mse(i) = mean(e.^2);
        disp("iter: [" + i + "], mse: " + mse(i) )
    end
    
    u = filter(w, 1, ref);

end

