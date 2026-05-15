function [u, w, y, mse] = fxrls(ref, target, ir, n_iter, filter_order, lamda, delta)

    w = zeros(filter_order, 1);
    y = zeros(1, length(target)); % the final output of (x_ref * w)[n]
    e = zeros(1, length(target)); % the error 
    
    x_filt = filter(ir, 1, ref);
    
    mse = zeros(n_iter, 1);

    % RLS state: inverse correlation matrix
    P = delta * eye(filter_order);
    
    for i = 1:n_iter
    
        x = zeros(filter_order, 1);
    
        for n = 1:length(target)
            x(2:end) = x(1:end-1);
            x(1) = x_filt(n);

            y(n) = w.' * x;
            e(n) = target(n) - y(n);
           
            % RLS update
            Px = P * x;
            denom = lamda + x.' * Px;     
            k = Px / denom; % gain vector (N x 1)

            w = w + k * e(n);
            P = (P - k * (x.' * P)) / lamda; % covariance update

        end
    
        mse(i) = mean(e.^2);
        disp("iter: [" + i + "], mse: " + mse(i) )
    end
    
    u = filter(w, 1, ref);

end

