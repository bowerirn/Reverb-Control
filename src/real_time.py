import numpy as np
import sounddevice as sd
from scipy.signal import lfilter

class RealtimeFxLMS:
    def __init__(
        self,
        audio_device,
        source,
        ir,
        block_size=1024,
        filter_order=256,
        eps=1e-8,
    ):
        self.ad = audio_device
        self.fs = self.ad.fs
        self.block_size = block_size
        self.M = filter_order
        self.eps = eps

        self.source = source
        if self.source.ndim > 1:
            self.source = self.source[:, 0]

        self.ir = ir

        self.reset()

    def reset(self):
        # Adaptive filter weights
        self.w = np.zeros(self.M, dtype=np.float32)

        # Delay lines
        self.x = np.zeros(self.M, dtype=np.float32)       # raw reference buffer
        self.xf = np.zeros(self.M, dtype=np.float32)      # filtered-reference buffer

        self.norm_sq = np.dot(self.xf, self.xf)

        # State for filtering ref through secondary path
        self.zi = np.zeros(len(self.ir) - 1, dtype=np.float32)

        N = len(self.source)
        self.error_log = np.zeros(N, dtype=np.float32)
        self.cancel_log = np.zeros(N, dtype=np.float32)
        self.source_log = np.zeros(N, dtype=np.float32)
        self.log_pos = 0
        
        num_blocks = int(np.ceil(len(self.source) / self.block_size))
        self.w_norm_log = np.zeros(num_blocks, dtype=np.float32)
        self.block_pos = 0

        self.pos = 0
        self.step = 0
        self.running = True


    def get_source_block(self, frames):
        end = self.pos + frames

        if end <= len(self.source):
            block = self.source[self.pos:end]
        else:
            block = np.zeros(frames, dtype=np.float32)
            remaining = max(0, len(self.source) - self.pos)
            if remaining > 0:
                block[:remaining] = self.source[self.pos:]

            self.running = False

        self.pos = end
        return block.astype(np.float32)



    def fxlms(
        self, 
        step_fn, 
        clean_source=True, 
        fx=True, 
        nlms=True, 
        leak=1e-5, 
        source_gain=1.0, 
        cancel_gain=0.3,
        max_norm=1.0
    ):
        
        def _nlms_update(
            ref_block, 
            cancel_block, 
            xf_block, 
            error_block, 
            frames, 
            step_size, 
            decay
        ):
            for n in range(frames):
                # Update raw reference delay line
                self.x[1:] = self.x[:-1] # is ring buffer faster here?
                self.x[0] = ref_block[n]


                # Controller output: y = w^T x
                y = np.dot(self.w, self.x)

                # Send opposite sign to cancel
                cancel_block[n] = -cancel_gain * y

                # Update filtered-reference delay line
                old = self.xf[-1]
                self.xf[1:] = self.xf[:-1]
                self.xf[0] = xf_block[n]

                e = error_block[n]

                self.norm_sq += self.xf[0]**2 - old**2
                norm = self.eps + self.norm_sq

                mu = step_size * e / norm
                self.w *= decay
                self.w += mu * self.xf

        def _lms_update(
            ref_block, 
            cancel_block, 
            xf_block, 
            error_block, 
            frames, 
            step_size, 
            decay
        ):
            for n in range(frames):
                self.x[1:] = self.x[:-1] # is ring buffer faster here?
                self.x[0] = ref_block[n]

                y = np.dot(self.w, self.x)

                cancel_block[n] = -cancel_gain * y

                self.xf[1:] = self.xf[:-1]
                self.xf[0] = xf_block[n]

                e = error_block[n]

                mu = step_size * e
                self.w *= decay
                self.w += mu * self.xf

        def callback(indata, outdata, frames, time, status):
            if status:
                print(status)

            source_block = self.get_source_block(frames) * source_gain

            error_block = indata[:, 0].astype(np.float32)

            ref_block = source_block if clean_source else indata[:, 1].astype(np.float32)

            if fx:
                xf_block, self.zi = lfilter(
                    self.ir,
                    np.array([1.0], dtype=np.float32),
                    ref_block,
                    zi=self.zi
                )
            else:
                xf_block = ref_block

            cancel_block = np.zeros(frames, dtype=np.float32)

            self.step += 1
            step_size = step_fn(self.step)

            decay = 1 - leak

            if nlms:
                _nlms_update(
                    ref_block, 
                    cancel_block, 
                    xf_block, 
                    error_block, 
                    frames, 
                    step_size, 
                    decay
                )
            else:
                _lms_update(
                    ref_block, 
                    cancel_block, 
                    xf_block, 
                    error_block, 
                    frames, 
                    step_size, 
                    decay
                )

            w_norm = np.linalg.norm(self.w)
            # if w_norm > max_norm:
            #     self.w *= max_norm / w_norm

            cancel_block = np.clip(cancel_block, -0.1, 0.1)

            end = min(self.log_pos + frames, len(self.error_log))
            ncopy = end - self.log_pos

            self.error_log[self.log_pos:end] = error_block[:ncopy]
            self.cancel_log[self.log_pos:end] = cancel_block[:ncopy]
            self.source_log[self.log_pos:end] = source_block[:ncopy]
            self.log_pos = end

            self.w_norm_log[self.block_pos] = w_norm
            self.block_pos += 1

            outdata[:, 0] = cancel_block
            outdata[:, 1] = source_block

        return callback

        

    def run(
        self, 
        clean_source=True, 
        fx=True, 
        nlms=True, 
        step_fn=0.01, 
        leak=1e-5, 
        source_gain=1.0, 
        cancel_gain=0.3,
        max_norm=1.0, 
    ):
        if isinstance(step_fn, (int, float)):
            step_size = float(step_fn)
            step_fn = lambda _: step_size

        with sd.Stream(
            samplerate=self.fs,
            blocksize=self.block_size,
            device=self.ad.device,
            channels=(2, 2),
            dtype="float32",
            callback=self.fxlms(
                step_fn, 
                clean_source=clean_source, 
                fx=fx, 
                nlms=nlms, 
                leak=leak, 
                source_gain=source_gain, 
                cancel_gain=cancel_gain,
                max_norm=max_norm, 
            ),
        ):
            while self.running:
                sd.sleep(100)

        source = self.source_log[:self.log_pos]
        cancel = self.cancel_log[:self.log_pos]
        error = self.error_log[:self.log_pos]
        w_norm_log = self.w_norm_log[:self.block_pos]

        return source, cancel, error, w_norm_log