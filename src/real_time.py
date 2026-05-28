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

        # State for filtering ref through secondary path
        self.zi = np.zeros(len(self.ir) - 1, dtype=np.float32)

        self.error_log = []
        self.cancel_log = []
        self.source_log = []
        self.w_norm_log = []

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



    def fxlms(self, step_fn, leak=1e-5, max_norm=1.0, source_gain=1.0, cancel_gain=0.3):

        def callback(indata, outdata, frames, time, status):
            if status:
                print(status)

            source_block = self.get_source_block(frames) * source_gain

            error_block = indata[:, 0].astype(np.float32)
            ref_block = indata[:, 1].astype(np.float32)
            # ref_block = source_block

            # Filter reference through estimated secondary path
            # xf_block, self.zi = lfilter(
            #     self.ir,
            #     [1.0],
            #     ref_block,
            #     zi=self.zi
            # )
            # xf_block = xf_block.astype(np.float32)
            xf_block = ref_block

            cancel_block = np.zeros(frames, dtype=np.float32)

            self.step += 1
            step_size = step_fn(self.step)

            decay = 1 - leak

            for n in range(frames):
                # Update raw reference delay line
                self.x[1:] = self.x[:-1]
                self.x[0] = ref_block[n]

                # Controller output: y = w^T x
                y = np.dot(self.w, self.x)

                # Send opposite sign to cancel
                cancel_block[n] = -cancel_gain * y

                # Update filtered-reference delay line
                self.xf[1:] = self.xf[:-1]
                self.xf[0] = xf_block[n]

                e = error_block[n]

                norm = self.eps + np.dot(self.xf, self.xf)
                self.w = self.w * decay + step_size / norm * e * self.xf

                # w_norm = np.linalg.norm(self.w)
                # if w_norm > max_norm:
                #     self.w *= max_norm / w_norm

            cancel_block = np.clip(cancel_block, -0.1, 0.1)

            self.error_log.append(error_block.copy())
            self.cancel_log.append(cancel_block.copy())
            self.source_log.append(source_block.copy())
            self.w_norm_log.append(np.linalg.norm(self.w))

            outdata[:, 0] = cancel_block
            outdata[:, 1] = source_block

        return callback

        

    def run(self, step_fn=0.01, leak=1e-5, max_norm=1.0, source_gain=1.0, cancel_gain=0.3):
        if isinstance(step_fn, (int, float)):
            step_size = float(step_fn)
            step_fn = lambda _: step_size

        with sd.Stream(
            samplerate=self.fs,
            blocksize=self.block_size,
            device=self.ad.device,
            channels=(2, 2),
            dtype="float32",
            callback=self.fxlms(step_fn, leak=leak, max_norm=max_norm, source_gain=source_gain, cancel_gain=cancel_gain),
        ):
            while self.running:
                sd.sleep(100)

        source = np.concatenate(self.source_log)
        cancel = np.concatenate(self.cancel_log)
        error = np.concatenate(self.error_log)

        return source, cancel, error