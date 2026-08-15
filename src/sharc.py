import matplotlib.pyplot as plt
from src.measure_ir import *
from pathlib import Path
import numpy as np
from src.audio_device import AudioDevice
from src.midi_protocol import MidiProtocol
from itertools import product


class Sharc:
    def __init__(self, source, ad: AudioDevice, midi_in=None, midi_out=None):
        self.ad = ad
        self.source = source
        self.no_cancels = {}
        self.midi_protocol = MidiProtocol(midi_in, midi_out)

        self.n_repeats = 1
        self.error_mic = []
        self.w_norm_log = []
        self.true_norm_log = []
        self.weights = np.array([])

        self.ir_file = r"..\experiments\data\current_irs.npz"
        self.panel_ir_h = r'..\sharc\reverb\reverb_core1\src'

        self.mu = None
        self.leak = None
        self.eps = None
        self.cancel_gain = None
        self.ref_threshold = None
        self.mavg_tau_ms = None
        self.lag = None
        self.update_sign = None

        if self.midi_protocol.connected:
            self.set(
                mu=1e-3,
                leak=3e-7,
                eps=1e-6,
                cancel_gain=0.1,
                ref_threshold=3e-4,
                mavg_tau_ms=100,
                lag=86,
                update_sign=1,
            )

    def set(
        self,
        mu=None,
        leak=None,
        eps=None,
        cancel_gain=None,
        ref_threshold=None,
        mavg_tau_ms=None,
        lag=None,
        update_sign=None,
    ):
        if mu is not None and mu != self.mu:
            self.set_mu(mu)
        if leak is not None and leak != self.leak:
            self.set_leak(leak)
        if eps is not None and eps != self.eps:
            self.set_eps(eps)
        if cancel_gain is not None and cancel_gain != self.cancel_gain:
            self.set_cancel_gain(cancel_gain)
        if ref_threshold is not None and ref_threshold != self.ref_threshold:
            self.set_ref_threshold(ref_threshold)
        if mavg_tau_ms is not None and mavg_tau_ms != self.mavg_tau_ms:
            self.set_mavg_tau_ms(mavg_tau_ms)
        if lag is not None and lag != self.lag:
            self.set_lag(lag)
        if update_sign is not None and update_sign != self.update_sign:
            self.set_update_sign(update_sign) 

    def close(self):
        self.midi_protocol.close()

    def set_mu(self, mu: float):
        self.midi_protocol.set_mu(mu)
        self.mu = mu

    def set_leak(self, leak: float):
        self.midi_protocol.set_leak(leak)
        self.leak = leak

    def set_eps(self, eps: float):
        self.midi_protocol.set_eps(eps)
        self.eps = eps

    def set_cancel_gain(self, gain: float):
        self.midi_protocol.set_cancel_gain(gain)
        self.cancel_gain = gain

    def set_ref_threshold(self, threshold: float):
        self.midi_protocol.set_ref_threshold(threshold)
        self.ref_threshold = threshold

    def set_mavg_tau_ms(self, tau_ms: int):
        self.midi_protocol.set_mavg_tau_ms(tau_ms)
        self.mavg_tau_ms = tau_ms

    def set_lag(self, lag: int):
        self.midi_protocol.set_lag(lag)
        self.lag = lag

    def set_update_sign(self, sign: int) -> None:
        self.midi_protocol.set_update_sign(sign)
        self.update_sign = sign

    def set_adapt(self, adapt: bool) -> None:
        self.midi_protocol.set_adapt(adapt)

    def set_off(self, off: bool) -> None:   
        self.midi_protocol.set_off(off)

    def seed_delta(self, index: int, amplitude: float = 1.0):
        self.midi_protocol.seed_delta(index, amplitude)

    def reset(self) -> None:
        self.midi_protocol.request_reset()
        self.w_norm_log = []
        self.error_log = []

    def get_weights(self) -> np.ndarray:
        weights, true_norm = self.midi_protocol.request_weights()
        self.weights = weights
        return weights

    def save_irs(self):
        error_ir, ref_ir = measure_ir(self.ad)
        np.savez(self.ir_file, error_ir=error_ir, ref_ir=ref_ir)

    def prep_irs(self, ir_len=256, panel_to_err_cm=4.5):

        irs = np.load(self.ir_file)
        error_ir = irs['error_ir']
        ref_ir = irs['ref_ir']
        
        error_ir, ref_ir = align_irs_by_distance(
            error_ir, ref_ir, 
            distance_cm=panel_to_err_cm, 
            ir_len=ir_len, 
            fs=self.ad.fs
        )

        plt.plot(error_ir, label='panel_ir')
        plt.plot(ref_ir, label='ref_ir')
        plt.legend()

        self.write_ir(error_ir, ref_ir, self.panel_ir_h)



    def cancel(
        self, n_repeats, adapt=True,
        mu=None,
        leak=None,
        eps=None,
        cancel_gain=None,
        ref_threshold=None,
        mavg_tau_ms=None,
        lag=None,
        update_sign=None,
    ):

        self.set(
            mu=mu,
            leak=leak,
            eps=eps,
            cancel_gain=cancel_gain,
            ref_threshold=ref_threshold,
            mavg_tau_ms=mavg_tau_ms,
            lag=lag,
            update_sign=update_sign,
        )
        self.set_off(False)
        self.set_adapt(adapt)

        source = np.tile(self.source, n_repeats)
        self.error_log, self.ref_log = self.ad.play(left=source)
        # for _ in range(n_repeats):
            # self.error_log, self.ref_log = self.ad.play(left=self.source)

            # weights, true_norm = self.midi_protocol.request_weights()
            # self.weights = weights
            # norm = np.linalg.norm(weights)

            # self.w_norm_log.append(norm)
            # self.true_norm_log.append(true_norm)

            # if (abs(true_norm - norm) > 1e-4):
                # print(f"Warning: weight norm is {true_norm:.6f}, which may indicate instability.")

        self.set_adapt(False)

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

            return db
        else:
            print(f"No baseline recording found for n_repeats={n_repeats}.")


    def grid_search(self, n_repeats, param_dict):
        # Scalars become one-element lists.
        search_values = {
            key: value if isinstance(value, (list, tuple, np.ndarray)) else [value]
            for key, value in param_dict.items()
        }

        keys = list(search_values.keys())

        best_params = None
        best_db = np.inf
        best_weights = None

        for combination in product(*(search_values[key] for key in keys)):
            params = dict(zip(keys, combination))

            # Only set parameters explicitly supplied to grid_search().
            self.set(**params)

            self.reset()
            print(f"Params: {params}")
            db = self.cancel(n_repeats)
            print()

            assert db is not None, f"No baseline available for n_repeats={n_repeats}"

            if db < best_db:
                best_db = db
                best_params = params.copy()
                # best_weights = self.weights.copy()

        print(
            f"Best Params: {best_params}, "
            f"Best ANC / No ANC: {best_db:.2f} dB"
        )

        return best_params, best_db#, best_weights


    def no_cancel(self, n_repeats):
        source = np.tile(self.source, n_repeats)

        self.set_adapt(False)
        self.set_off(True)
        error_mic, self.ref_log = self.ad.play(left=source)
        self.set_off(False)
        
        self.no_cancels[n_repeats] = error_mic

    
    def write_ir(self, error_ir, ref_ir, path='.'):
        path = Path(path)

        if path.is_dir() or path.suffix == "":
            path = path / "panel_ir.h"

        path.parent.mkdir(parents=True, exist_ok=True)

        with path.open("w") as f:
            f.write("#pragma once\n\n")
            f.write(f"#define IR_LEN {len(error_ir)}\n\n")

            f.write("static const float error_ir[IR_LEN] = {\n")
            for v in error_ir:
                f.write(f"    {float(v):.9e}f,\n")
            f.write("};\n\n")

            f.write("static const float ref_ir[IR_LEN] = {\n")
            for v in ref_ir:
                f.write(f"    {float(v):.9e}f,\n")
            f.write("};\n")




    def plot_weights(self, weights=None, title_ext=''):
        if weights is None:
            weights = self.weights

        plt.figure(figsize=(10, 4))
        plt.plot(weights)
        plt.title(f"Adaptive Filter Weights {title_ext}")
        plt.xlabel("Weight Index")
        plt.ylabel("Weight Value")
        plt.grid(True)
        plt.show()


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
        self.plot_w_norm(title_ext)
