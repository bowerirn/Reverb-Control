import numpy as np
import sounddevice as sd
from scipy.signal import lfilter

class RealtimeFxLMS:
    def __init__(
        self,
        audio_device,
        source,
        ir,
        blocksize=1024,
        filter_order=256,
        step_fn=0.01,
        eps=1e-8,
        source_gain=0.5,
        cancel_gain=1.0,
    ):
        self.ad = audio_device
        self.fs = self.ad.fs
        self.blocksize = blocksize
        self.M = filter_order
        self.eps = eps
        self.source_gain = source_gain
        self.cancel_gain = cancel_gain

        if isinstance(step_fn, (int, float)):
            self.step_fn = lambda i: step_fn
        else: self.step_fn = step_fn

        self.source = np.asarray(source, dtype=np.float32).squeeze()
        if self.source.ndim > 1:
            self.source = self.source[:, 0]

        self.ir = np.asarray(ir, dtype=np.float32).squeeze()

        # Adaptive filter weights
        self.w = np.zeros(self.M, dtype=np.float32)

        # Delay lines
        self.x = np.zeros(self.M, dtype=np.float32)       # raw reference buffer
        self.xf = np.zeros(self.M, dtype=np.float32)      # filtered-reference buffer

        # State for filtering ref through secondary path
        self.zi = np.zeros(len(self.ir) - 1, dtype=np.float32)

        self.pos = 0
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

    def callback(self, indata, outdata, frames, time, status):
        if status:
            print(status)

        source_block = self.get_source_block(frames) * self.source_gain

        error_block = indata[:, 0].astype(np.float32)
        ref_mic_block = indata[:, 1].astype(np.float32)

        # Filter reference through estimated secondary path
        xf_block, self.zi = lfilter(
            self.ir,
            [1.0],
            ref_mic_block,
            zi=self.zi
        )
        xf_block = xf_block.astype(np.float32)

        cancel_block = np.zeros(frames, dtype=np.float32)

        for n in range(frames):
            # Update raw reference delay line
            self.x[1:] = self.x[:-1]
            self.x[0] = ref_mic_block[n]

            # Controller output: y = w^T x
            y = np.dot(self.w, self.x)

            # Send opposite sign to cancel
            cancel_block[n] = -self.cancel_gain * y

            # Update filtered-reference delay line
            self.xf[1:] = self.xf[:-1]
            self.xf[0] = xf_block[n]

            e = error_block[n]

            norm = self.eps + np.dot(self.xf, self.xf)
            self.w += (self.step_fn(n) / norm) * e * self.xf

        outdata[:, 0] = cancel_block
        outdata[:, 1] = source_block

    def run(self):
        with sd.Stream(
            samplerate=self.fs,
            blocksize=self.blocksize,
            device=self.ad.device,
            channels=(2, 2),
            dtype="float32",
            callback=self.callback,
        ):
            while self.running:
                sd.sleep(100)