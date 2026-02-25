classdef FD_FxLMS < AdaptiveFilter
    properties
        eps = 1e-8;  % small constant to avoid division by zero
        step_size = 0.05;
        block_len = 128; % number of samples per block (B)

        W;
        Xhist;
        Pow;
        alpha;  
        x_ov;
        y_ov;
        xpad;
        tpad;
        nBlocks;
        padN;
        ypad;
        epad;
    end

    methods
        function obj = FD_FxLMS(ref, target, ir, filter_order, step_size, block_len, varargin)
            obj@AdaptiveFilter(ref, target, ir, filter_order, 'FD_FxLMS');
            
            obj.step_size = step_size;
            obj.block_len = block_len;
            if nargin > 6
                obj.eps = varargin{1};
            end

            % Partitioning
            P = ceil(filter_order / block_len);          % # partitions
            Nfft = 2 * block_len;             % common choice
            if Nfft < 2*block_len
                error("Nfft must be >= 2*block_len");
            end

            % Frequency-domain adaptive filter partitions Wp(k)
            obj.W = complex(zeros(Nfft, P));   % each column is one partition in freq-domain

            % Keep last (P-1) FFT blocks of input spectra for convolution
            obj.Xhist = complex(zeros(Nfft, P));  % newest in col 1

            % Power estimate per bin for normalization (smoothed)
            obj.Pow = eps * ones(Nfft, 1);
            obj.alpha = 0.9; % smoothing for power estimate (0.8-0.99 typical)

            % Overlap-save buffers
            obj.x_ov = zeros(block_len, 1);   % previous block tail (for building 2*block_len input)
            obj.y_ov = zeros(block_len, 1);   % overlap part from ifft result (discarded)

            % Pad signals to a multiple of B
            obj.nBlocks = ceil(length(target) / block_len);
            obj.padN = obj.nBlocks * block_len - length(target);
            obj.xpad = [x_filt; zeros(padN, 1)];
            obj.tpad = [target; zeros(padN, 1)];
        end

        
        function to_gpu(obj); end
        function to_cpu(obj); end

        
        function update_coefficients(obj)
            obj.Xhist(:) = 0;
            obj.x_ov(:) = 0;


            for b = 1:obj.nBlocks
                idx = (b-1)*obj.block_len + (1:obj.block_len);

                % Build 2B block for FFT (overlap-save)
                x_blk = obj.xpad(idx);
                x2 = [obj.x_ov; x_blk];         % length 2B
                obj.x_ov = x_blk;               % next overlap

                X = fft(x2, Nfft);

                % Update input spectrum history (shift right, insert newest)
                obj.Xhist(:,2:end) = obj.Xhist(:,1:end-1);
                obj.Xhist(:,1) = X;

                % Compute output spectrum: Y = sum_p Wp .* Xhist_p
                Y = zeros(Nfft, 1);
                for p = 1:P
                    Y = Y + obj.W(:,p) .* obj.Xhist(:,p);
                end

                y2 = real(ifft(Y, Nfft));   % time-domain 2B (valid part is last B)
                y_blk = y2(obj.block_len+1 : obj.block_len+obj.block_len);

                % Error block (time domain)
                d_blk = obj.tpad(idx);
                e_blk = d_blk - y_blk;

                % Store
                obj.ypad(idx) = y_blk;
                obj.epad(idx) = e_blk;

                % FFT of error, zero-pad as 2B with leading zeros (alignment)
                E = fft([zeros(B,1); e_blk], Nfft);

                % Power estimate per bin for normalization (smoothed)
                obj.Pow = obj.alpha * obj.Pow + (1 - obj.alpha) * (abs(X).^2);
                
                % Block-scaled step size (FD updates are "chunkier")
                mu_fd = obj.step_size / obj.block_len;
                
                % Safer normalization (avoid tiny denominators)
                den = obj.eps + obj.Pow;
                den = max(den, 1e-8);
                
                % Gradient term
                G = (mu_fd ./ den) .* E;  % Nfft x 1
                
                for p = 1:P
                    % Update in frequency domain
                    obj.W(:,p) = obj.W(:,p) + conj(obj.Xhist(:,p)) .* G;
                
                    % ---- IMPORTANT: enforce time-domain partition constraint ----
                    wp = real(ifft(obj.W(:,p), Nfft));
                    wp(B+1:end) = 0;      % keep only first B taps of this partition
                    obj.W(:,p) = fft(wp, Nfft);
                end

            end

        % Trim padding
        obj.y = obj.ypad(1:N);
        obj.e = obj.epad(1:N);
        end

    end
end