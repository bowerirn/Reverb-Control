import matplotlib.pyplot as plt
from src.measure_ir import *
from pathlib import Path
import numpy as np


class Sharc:
    def __init__(self, source, ad: AudioDevice):
        self.ad = ad
        self.source = source
        self.no_cancels = {}


        self.n_repeats = 1
        self.error_mic = []



    def prep_ir(self, ir_len=128, panel_to_err_cm=7):
        error_ir, _ = measure_ir(self.ad)
        panel_ir = align_ir_by_distance(error_ir, panel_to_err_cm, ir_len)
        plt.plot(panel_ir)
        self.write_ir(panel_ir, r'..\sharc\reverb\reverb_core1\src')



    def cancel(self, n_repeats):
        source = np.tile(self.source, n_repeats)

        self.error_log, self.ref_log = self.ad.play(left=source)
        self.n_repeats = n_repeats

        if n_repeats in self.no_cancels:
            err = self.error_log
            nc = self.no_cancels[n_repeats]

            n = min(len(err), len(nc))
            err = err[:n]
            nc = nc[:n]

            db = 10 * np.log10(
                (np.mean(err**2) + 1e-20) /
                (np.mean(nc**2) + 1e-20)
            )

            print(f"Overall ANC / No ANC: {db:.2f} dB")
        else:
            print(f"No baseline recording found for n_repeats={n_repeats}.")


    def no_cancel(self, n_repeats):
        source = np.tile(self.source, n_repeats)
        error_mic, self.ref_log = self.ad.play(left=source)
        self.no_cancels[n_repeats] = error_mic

    
    def write_ir(self, ir, path='.'):
        path = Path(path)

        if path.is_dir() or path.suffix == "":
            path = path / "panel_ir.h"

        path.parent.mkdir(parents=True, exist_ok=True)

        with path.open("w") as f:
            f.write("#pragma once\n\n")
            f.write(f"#define IR_LEN {len(ir)}\n\n")
            f.write("static const float panel_ir[IR_LEN] = {\n")
            for v in ir:
                f.write(f"    {float(v):.9e}f,\n")
            f.write("};\n")





    def plot_error_mic(self, title_ext=''):
        plt.figure(figsize=(10, 4))

        nc = self.no_cancels.get(self.n_repeats, None)
        if nc is not None:
            plt.plot(nc, alpha=0.8, label="No Cancel")

        plt.plot(self.error_log, alpha=0.8, label="With Cancel")


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
    
    def plot_rms(self, title_ext=''):
        nc = self.no_cancels.get(self.n_repeats, None)
        if nc is not None:
            err_rms_nc = self._block_rms(nc)

        err_rms = self._block_rms(self.error_log)

        plt.figure(figsize=(8, 4))

        if nc is not None:
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



    def plot_cumulative_error_reduction(self, skip_s=1.0):
        err = self.error_log
        nc = self.no_cancels.get(self.n_repeats, None)

        n = min(len(err), len(nc))
        err = err[:n]
        nc = nc[:n]

        skip = int(skip_s * self.ad.fs)
        err = err[skip:]
        nc = nc[skip:]

        e2 = np.cumsum(err**2)
        n2 = np.cumsum(nc**2)

        reduction_db = 10 * np.log10((e2 + 1e-20) / (n2 + 1e-20))
        t = np.arange(len(err)) / self.ad.fs + skip_s

        plt.figure(figsize=(8, 4))
        plt.plot(t, reduction_db)
        plt.axhline(0, ls="--")
        plt.grid(True)

        plt.xlabel("Time [s]")
        plt.ylabel("Cumulative ANC / No ANC [dB]")
        plt.title("Cumulative Error Reduction Over Time")
        plt.show()



    def plot_error_ratio(self, skip_s=1.0, window_s=0.25):
        err = np.asarray(self.error_log)
        nc = np.asarray(self.no_cancels.get(self.n_repeats))

        n = min(len(err), len(nc))
        err = err[:n]
        nc = nc[:n]

        fs = self.ad.fs
        skip = int(skip_s * fs)

        err = err[skip:]
        nc = nc[skip:]

        window = max(1, int(window_s * fs))
        kernel = np.ones(window) / window

        # Moving mean-square energy
        err_power = np.convolve(err**2, kernel, mode="valid")
        nc_power = np.convolve(nc**2, kernel, mode="valid")

        relative_db = 10 * np.log10(
            (err_power + 1e-20) / (nc_power + 1e-20)
        )

        # Time at the center of each window
        t = (
            np.arange(len(relative_db))
            + window / 2
        ) / fs + skip_s

        plt.figure(figsize=(8, 4))
        plt.plot(t, relative_db)
        plt.axhline(0, linestyle="--")
        plt.grid(True)

        plt.xlabel("Time [s]")
        plt.ylabel("ANC / No ANC [dB]")
        plt.title(f"Error Level Relative to No ANC ({window_s:.2f} s window)")
        plt.show()


    def all_plots(self, title_ext=''):
        self.plot_error_mic(title_ext)

        nc = self.no_cancels.get(self.n_repeats, None)
        if nc is not None:
            self.plot_error_reduction()
            self.plot_cumulative_error_reduction()
        # self.plot_w_norm(title_ext)
