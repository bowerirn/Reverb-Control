import numpy as np
import sounddevice as sd
from scipy.signal import lfilter, coherence, welch, correlate
import matplotlib.pyplot as plt

from . import lms_rt_ext

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

        self.delay_ms = None
        self.reset()

    def reset(self, lag=0):
        self.adapt_delay = lag
        self.control_delay = self.adapt_delay   # TODO: calculate true delay here?

        # Adaptive filter weights
        self.w = np.zeros(self.M, dtype=np.float32)

        # Delay lines
        self.x = np.zeros(self.M + self.control_delay, dtype=np.float32)       # raw reference buffer
        self.xf = np.zeros(self.M + self.adapt_delay, dtype=np.float32)      # filtered-reference buffer

        self.norm_sq = np.dot(self.xf, self.xf)

        # State for filtering ref through secondary path
        self.zi = np.zeros(len(self.ir) - 1, dtype=np.float32)

        self.heads = np.zeros(2, dtype=np.int32) # for lms ring buffers in cpp

        N = len(self.source)
        self.error_log = np.zeros(N, dtype=np.float32)
        self.cancel_log = np.zeros(N, dtype=np.float32)
        self.log_pos = 0
        
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



    def fxlms(
        self, 
        step_fn, 
        clean_source=True, 
        fx=True, 
        nlms=True, 
        leak=1e-5, 
        source_gain=1.0, 
        cancel_gain=0.3,
        w_update_sign=1.0,
        max_norm=1.0,
    ):
        
        
        
        # def _nlms_update(
        #     ref_block, 
        #     cancel_block, 
        #     xf_block, 
        #     error_block, 
        #     frames, 
        #     step_size, 
        #     decay
        # ):
        #     for n in range(frames):
        #         # Update raw reference delay line
        #         self.x[1:] = self.x[:-1] # is ring buffer faster here?
        #         self.x[0] = ref_block[n]


        #         # Controller output: y = w^T x
        #         y = np.dot(self.w, self.x)

        #         # Send opposite sign to cancel
        #         cancel_block[n] = -cancel_gain * y

        #         # Update filtered-reference delay line
        #         old = self.xf[-1]
        #         self.xf[1:] = self.xf[:-1]
        #         self.xf[0] = xf_block[n]

        #         e = error_block[n]

        #         self.norm_sq += self.xf[0]**2 - old**2
        #         norm = self.eps + self.norm_sq

        #         mu = step_size * e / norm
        #         self.w *= decay
        #         self.w += mu * self.xf

        # def _lms_update(
        #     ref_block, 
        #     cancel_block, 
        #     xf_block, 
        #     error_block, 
        #     frames, 
        #     step_size, 
        #     decay
        # ):
        #     for n in range(frames):
        #         self.x[1:] = self.x[:-1] # is ring buffer faster here?
        #         self.x[0] = ref_block[n]

        #         y = np.dot(self.w, self.x)

        #         cancel_block[n] = -cancel_gain * y

        #         self.xf[1:] = self.xf[:-1]
        #         self.xf[0] = xf_block[n]

        #         e = error_block[n]

        #         mu = step_size * e
        #         self.w *= decay
        #         self.w += mu * self.xf

        if max_norm is None:
            max_norm = 0.0

        def callback(indata, outdata, frames, time, status):
            self.delay_ms = time.outputBufferDacTime - time.inputBufferAdcTime
            

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

            update = lms_rt_ext.nlms_realtime_update if nlms else lms_rt_ext.lms_realtime_update
            w_norm = update(
                ref_block,
                xf_block.astype(np.float32, copy=False),
                error_block,
                self.x,
                self.xf,
                cancel_block,
                self.w,
                self.heads,
                float(step_size),
                float(leak),
                float(cancel_gain),
                float(w_update_sign),
                float(self.eps),
                float(max_norm),
                int(self.adapt_delay),
                int(self.control_delay),
            )

            cancel_block = np.clip(cancel_block, -0.1, 0.1)

            end = min(self.log_pos + frames, len(self.error_log))
            ncopy = end - self.log_pos

            self.error_log[self.log_pos:end] = error_block[:ncopy]
            self.cancel_log[self.log_pos:end] = cancel_block[:ncopy]
            self.log_pos = end

            self.w_norm_log.append(w_norm)

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
        w_update_sign=1.0,
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
                w_update_sign=w_update_sign,
                max_norm=max_norm, 
            ),
        ):
            while self.running:
                sd.sleep(100)

        cancel = self.cancel_log[:self.log_pos]
        error = self.error_log[:self.log_pos]
        w_norm_log = np.asarray(self.w_norm_log)

        title_ext = f'{"{"}step={step_fn(0):1.1e}, c_gain={cancel_gain}, leak={leak:1.1e}{", fx" if fx else ""}{", nlms" if nlms else ""}{", source" if clean_source else ''}{"}"}'

        return error, cancel, w_norm_log, title_ext
    

    def empirical_delay(self):
        err = self.error_log[:self.log_pos]
        src = self.source[:self.log_pos]

        corr = correlate(err, src, mode='full')
        lags = np.arange(-len(src)+1, len(err))

        lag = lags[np.argmax(np.abs(corr))]
        print("Lag from cross-correlation:", lag)
        print('Delay from device:', round(self.delay_ms * self.ad.fs))
        return lag
    
    

    def plot_error_mic(self, error_nc=None, title_ext=''):
        plt.figure(figsize=(10, 4))

        if error_nc is not None:
            plt.plot(error_nc, alpha=0.8, label="No Cancel")

        plt.plot(self.error_log[:self.log_pos], alpha=0.8, label="With Cancel")


        plt.title(f"Error Mic Signal During FxLMS {title_ext}")
        plt.xlabel("Samples")
        plt.ylabel("Amplitude")
        plt.legend(loc='upper right')
        plt.grid(True)
        plt.show()


    def _block_rms(self, x):
        n = len(x) // self.block_size
        x = x[:n * self.block_size].reshape(n, self.block_size)
        return np.sqrt(np.mean(x**2, axis=1))
    
    def plot_rms(self, error_nc=None, title_ext=''):
        if error_nc is not None:
            err_rms_nc = self._block_rms(error_nc)

        err_rms = self._block_rms(self.error_log[:self.log_pos])

        plt.figure(figsize=(8, 4))

        if error_nc is not None:
            plt.plot(err_rms_nc, label="No Cancel")

        plt.plot(err_rms, label="With Cancel")

        plt.yscale("log")
        plt.title(f"Error Mic RMS Over Time {title_ext}")
        plt.xlabel("Block")
        plt.ylabel("RMS")
        plt.legend()
        plt.grid(True)
        plt.show()


    def plot_w_norm(self, title_ext=''):
        plt.figure(figsize=(8, 4))
        plt.plot(np.asarray(self.w_norm_log))
        plt.xlabel("Update Step")
        plt.ylabel(r"$||w||_2$")
        plt.title(f"Adaptive Filter Weight Norm {title_ext}")
        plt.grid(True)
        plt.show()


    def plot_coherence(self, delay=None):
        err = self.error_log[:self.log_pos]
        src = self.source[:self.log_pos]

        if delay is None:
            delay = self.ad.delay_samples
        
        if delay > 0:
            src = src[:-delay]
            err = err[delay:]

        n = min(len(src), len(err))
        src = src[:n]
        err = err[:n]

        f, Cxy = coherence(src, err, fs=self.ad.fs)

        plt.figure(figsize=(8,4))
        plt.semilogx(f, Cxy)
        plt.ylim([0, 1])
        plt.grid(True, which="both")
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Coherence")
        plt.title(f"Reference → Error Mic Coherence, delay={delay}")
        plt.show()


    def plot_error_reduction(self, error_nc):
        f, Pc = welch(self.error_log[:self.log_pos], self.ad.fs)
        _, Pnc = welch(error_nc, self.ad.fs)

        reduction_db = 10 * np.log10(
            (Pc + 1e-20) /
            (Pnc + 1e-20)
        )

        plt.figure(figsize=(8,4))
        plt.semilogx(f, reduction_db)

        plt.axhline(0, ls='--')
        plt.grid(True, which="both")

        plt.xlabel("Frequency [Hz]")
        plt.ylabel("ANC / No ANC [dB]")
        plt.title("Spectral Error Reduction")

        plt.show()


    def plot_psd(self, error_nc=None):
        f, Ps = welch(self.source, self.ad.fs)
        _, Pc = welch(self.error_log[:self.log_pos], self.ad.fs)

        plt.figure(figsize=(8,4))

        plt.semilogx(f, 10 * np.log10(Ps + 1e-20), label="Source")
        plt.semilogx(f, 10 * np.log10(Pc + 1e-20), label="Error ANC")

        if error_nc is not None:
            _, Pnc = welch(error_nc, self.ad.fs)
            plt.semilogx(
                f, 10 * np.log10(Pnc + 1e-20), label="Error No ANC")

        plt.grid(True, which="both")
        plt.legend()

        plt.xlabel("Frequency [Hz]")
        plt.ylabel("PSD [dB]")
        plt.title("Signal Spectra")

        plt.show()


    def all_plots(self, error_nc=None, title_ext=''):
        self.plot_error_mic(error_nc, title_ext)
        if error_nc is not None:
            self.plot_error_reduction(error_nc)
        self.plot_psd(error_nc)
        self.plot_w_norm(title_ext)


    