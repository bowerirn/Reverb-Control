function [u, w, y, mse] = fd_fxlms(ref, target, ir, n_iter, filter_order, mu, eps, block_len)
% Frequency-domain (partitioned-block) FXLMS
%
% ref, target: same length signals
% ir: secondary-path estimate (used to prefilter ref -> x_filt)
% n_iter: number of passes over the data (like your code)
% filter_order: adaptive FIR length (L)
% mu, eps: NLMS params
% block_len: block size (B). Typical: 128, 256, 512

    ref = ref(:);
    target = target(:);
    N = length(target);
    L = filter_order;
    B = block_len;

    % Filtered reference (Fx) (same as your time-domain version)
    x_filt = filter(ir, 1, ref);

    % Partitioning
    P = ceil(L / B);          % # partitions
    Nfft = 2 * B;             % common choice
    if Nfft < 2*B
        error("Nfft must be >= 2*B");
    end

    % Frequency-domain adaptive filter partitions Wp(k)
    W = complex(zeros(Nfft, P));   % each column is one partition in freq-domain

    % Keep last (P-1) FFT blocks of input spectra for convolution
    Xhist = complex(zeros(Nfft, P));  % newest in col 1

    % Power estimate per bin for normalization (smoothed)
    Pow = eps * ones(Nfft, 1);
    alpha = 0.9; % smoothing for power estimate (0.8-0.99 typical)

    y = zeros(N, 1);
    e = zeros(N, 1);
    mse = zeros(n_iter, 1);

    % Overlap-save buffers
    x_ov = zeros(B, 1);   % previous block tail (for building 2B input)
    y_ov = zeros(B, 1);   % overlap part from ifft result (discarded)

    % Pad signals to a multiple of B
    nBlocks = ceil(N / B);
    padN = nBlocks * B - N;
    xpad = [x_filt; zeros(padN, 1)];
    tpad = [target; zeros(padN, 1)];

    for it = 1:n_iter
        % reset overlap/history each epoch to match your "full pass" behavior
        Xhist(:) = 0;
        x_ov(:) = 0;

        ypad = zeros(nBlocks*B, 1);
        epad = zeros(nBlocks*B, 1);

        for b = 1:nBlocks
            idx = (b-1)*B + (1:B);

            % Build 2B block for FFT (overlap-save)
            x_blk = xpad(idx);
            x2 = [x_ov; x_blk];         % length 2B
            x_ov = x_blk;               % next overlap

            X = fft(x2, Nfft);

            % Update input spectrum history (shift right, insert newest)
            Xhist(:,2:end) = Xhist(:,1:end-1);
            Xhist(:,1) = X;

            % Compute output spectrum: Y = sum_p Wp .* Xhist_p
            Y = zeros(Nfft, 1);
            for p = 1:P
                Y = Y + W(:,p) .* Xhist(:,p);
            end

            y2 = real(ifft(Y, Nfft));   % time-domain 2B (valid part is last B)
            y_blk = y2(B+1 : B+B);

            % Error block (time domain)
            d_blk = tpad(idx);
            e_blk = d_blk - y_blk;

            % Store
            ypad(idx) = y_blk;
            epad(idx) = e_blk;

            % FFT of error, zero-pad as 2B with leading zeros (alignment)
            E = fft([zeros(B,1); e_blk], Nfft);

            % Power estimate per bin for normalization (smoothed)
            Pow = alpha * Pow + (1 - alpha) * (abs(X).^2);
            
            % Block-scaled step size (FD updates are "chunkier")
            mu_fd = mu / B;
            
            % Safer normalization (avoid tiny denominators)
            den = eps + Pow;
            den = max(den, 1e-8);
            
            % Gradient term
            G = (mu_fd ./ den) .* E;  % Nfft x 1
            
            for p = 1:P
                % Update in frequency domain
                W(:,p) = W(:,p) + conj(Xhist(:,p)) .* G;
            
                % ---- IMPORTANT: enforce time-domain partition constraint ----
                wp = real(ifft(W(:,p), Nfft));
                wp(B+1:end) = 0;      % keep only first B taps of this partition
                W(:,p) = fft(wp, Nfft);
            end

        end

        % Trim padding
        y = ypad(1:N);
        e = epad(1:N);

        mse(it) = mean(e.^2);
        disp("iter: [" + it + "], mse: " + mse(it))
    end

    % Convert partitioned W back to time-domain FIR taps w (length L)
    w_full = zeros(P*B, 1);
    for p = 1:P
        wp_time = real(ifft(W(:,p), Nfft));
        wp_time = wp_time(1:B);                % first B samples are the partition taps
        w_full((p-1)*B + (1:B)) = wp_time;
    end
    w = w_full(1:L);

    % Your original returns "u = filter(w,1,ref)" (control from raw ref)
    u = filter(w, 1, ref);
end
