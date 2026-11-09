import numpy as np
import matplotlib.pyplot as plt
import math


class OfflineFxLMS:
    def __init__(
        self,
        M,
        error_ir,
        ref_ir,
        fs=96_000,
        nlms=True,
        dtype=np.float32,
        nonlinear=False,
    ):
        self.nonlinear = bool(nonlinear)
        self.M = int(M)
        self.path_ir = np.asarray(error_ir, dtype=dtype)
        self.feedback_ir = np.asarray(ref_ir, dtype=dtype)
        self.IR_LENGTH = len(self.path_ir)
        self.NLMS = bool(nlms)
        self.dtype = dtype
        self.fs = fs

        if len(self.feedback_ir) != self.IR_LENGTH:
            raise ValueError("error_ir and ref_ir must have the same length")

        # Public parameters -- same defaults as SHARC
        self.mu = dtype(1e-4)
        self.eps = dtype(1e-6)
        self.leak = dtype(3e-7)
        self.cancel_gain = dtype(0.02)
        self.update_sign = dtype(1.0)

        self.lag = 86

        

        self.adapt = False

        self.ref_threshold = dtype(3e-4)
        self.mavg_tau_ms = dtype(100)
        self.mavg_weight = np.exp(-1.0 / (self.mavg_tau_ms * 0.001 * self.fs))

        self.reset()


    def reset(self):
        M = self.M
        L = self.IR_LENGTH
        dtype = self.dtype

        self.x_head = 0
        self.xf_head = 0
        self.path_ir_head = 0
        self.feedback_ir_head = 0

        self.xnorm = dtype(0.0)
        self.mavg = dtype(0.0)

        self.update = False

        # Diagnostics
        self.max_step = dtype(0.0)
        self.min_xnorm = dtype(1e30)
        self.max_control = dtype(0.0)

        # In nonlinear mode, rows correspond to x, x**2 and x**3.
        branches = (3,) if self.nonlinear else ()

        # Filter weights
        self.w = np.zeros(branches + (M,), dtype=dtype)

        # Duplicated ring buffers
        self.x = np.zeros(branches + (2 * M,), dtype=dtype)

        self.XF_LEN = M + self.lag
        self.xf = np.zeros(branches + (2 * self.XF_LEN,), dtype=dtype)


        self.z_path = np.zeros(branches + (2 * L,), dtype=dtype)
        self.z_feedback = np.zeros(2 * L, dtype=dtype)

        self.xnorm_history = np.zeros(M, dtype=self.dtype)


    def set(
        self,
        mu=None,
        eps=None,
        leak=None,
        cancel_gain=None,
        update_sign=None,
        lag=None,
        ref_threshold=None,
        mavg_tau_ms=None,
        nonlinear=None,
    ):
        if mu is not None:
            self.mu = self.dtype(mu)
        if eps is not None:
            self.eps = self.dtype(eps)
        if leak is not None:
            self.leak = self.dtype(leak)
        if cancel_gain is not None:
            self.cancel_gain = self.dtype(cancel_gain)
        if update_sign is not None:
            self.update_sign = self.dtype(update_sign)
        mode_changed = nonlinear is not None and bool(nonlinear) != self.nonlinear
        if mode_changed:
            self.nonlinear = bool(nonlinear)
        if lag is not None:
            self.lag = int(lag)
        if lag is not None or mode_changed:
            self.reset()
        if ref_threshold is not None:
            self.ref_threshold = self.dtype(ref_threshold)
        if mavg_tau_ms is not None:
            self.mavg_tau_ms = int(mavg_tau_ms)
            self.mavg_weight = np.exp(-1.0 / (self.mavg_tau_ms * 0.001 * self.fs))
    


    # ------------------------------------------------------------------
    # Duplicated-ring helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _ring_push(new_sample, state, head, N):
        head -= 1

        if head < 0:
            head = N - 1

        old = state[head] if state.ndim == 1 else state[:, head].copy()
        state[..., head] = new_sample
        state[..., head + N] = new_sample

        return old, head

    @staticmethod
    def _ring_dot(coeffs, state, head, N):
        contiguous_state = state[head:head + N]
        return np.dot(coeffs, contiguous_state)

    @classmethod
    def _ring_fir(cls, new_sample, coeffs, state, head, N):

        _, head = cls._ring_push(new_sample, state, head, N)
        y = cls._ring_dot(coeffs, state, head, N)

        return y, head




    def seed_filter(self, coeffs, negate=False):
        coeffs = np.asarray(coeffs, dtype=self.dtype)

        if self.nonlinear:
            self.w.fill(0)
            if coeffs.ndim == 2:
                if coeffs.shape[0] != 3:
                    raise ValueError("Expected three rows for linear, square and cube weights")
                n = min(coeffs.shape[1], self.M)
                self.w[:, :n] = -coeffs[:, :n] if negate else coeffs[:, :n]
                return
        weights = self.w[0] if self.nonlinear else self.w
        n = min(len(coeffs), self.M)

        weights[:n] = -coeffs[:n] if negate else coeffs[:n]
        weights[n:] = 0.0

    def seed_delta(self, delay, amp):
        if 0 <= delay < self.M:
            weights = self.w[0] if self.nonlinear else self.w
            weights[delay] = self.dtype(amp)




    def _features(self, sample):
        x = self.dtype(sample)
        square = x * x
        return (x, square, square * x)

    def compute_control(self, ref, clean_feedback=False):
        cleaned_ref = ref
        if clean_feedback:
            predicted_feedback = self._ring_dot(self.feedback_ir, self.z_feedback, self.feedback_ir_head, self.IR_LENGTH)
            cleaned_ref -= predicted_feedback

        if self.nonlinear:
            _, self.x_head = self._ring_push(self._features(cleaned_ref), self.x, self.x_head, self.M)
            history = self.x[:, self.x_head:self.x_head + self.M]
            control_raw = np.einsum("bk,bk->", self.w, history)
        else:
            control_raw, self.x_head = self._ring_fir(cleaned_ref, self.w, self.x, self.x_head, self.M)
        control = -self.cancel_gain * control_raw

        if clean_feedback:
            _, self.feedback_ir_head = self._ring_push(control, self.z_feedback, self.feedback_ir_head, self.IR_LENGTH)

        abs_control = abs(control)
        
        if abs_control > self.max_control:
            self.max_control = abs_control

        return control, cleaned_ref
    
    def process(self, cleaned_ref, error_mic):

        if not self.adapt or self.mu == 0.0:
            return

        if self.nonlinear:
            _, self.path_ir_head = self._ring_push(
                self._features(cleaned_ref), self.z_path, self.path_ir_head, self.IR_LENGTH
            )
            history = self.z_path[:, self.path_ir_head:self.path_ir_head + self.IR_LENGTH]
            xf_sample = history @ self.path_ir
        else:
            xf_sample, self.path_ir_head = self._ring_fir(cleaned_ref, self.path_ir, self.z_path, self.path_ir_head, self.IR_LENGTH)

        old, self.xf_head = self._ring_push(xf_sample, self.xf, self.xf_head, self.XF_LEN)

        
        update_head = self.xf_head + self.lag

        if update_head >= self.XF_LEN:
            update_head -= self.XF_LEN
        

        if self.nonlinear:
            newest = self.xf[:, update_head]
            self.xnorm += np.dot(newest, newest) - np.dot(old, old)
        else:
            newest = self.xf[update_head]
            self.xnorm += newest * newest - old * old

        if self.xnorm < 0.0:
            self.xnorm = 0.0


        self.mavg = self.mavg_weight * self.mavg + (1.0 - self.mavg_weight) * self.xnorm
        if cleaned_ref < .001:
            return       
        # if self.mavg < self.ref_threshold:
        #     return 



        step = self.mu

        if self.NLMS:
            step = self.mu / (self.eps + self.xnorm)


        self.update = not self.update
        if not self.update:
            return 

        
        update_scale = self.update_sign * step * error_mic
        decay = 1.0 - self.leak


        xf_contiguous = self.xf[..., update_head:update_head + self.M]


        if self.xnorm < self.min_xnorm:
            self.min_xnorm = self.xnorm

        if step > self.max_step:
            self.max_step = step

        self.debug_xf = xf_contiguous.copy()
        self.debug_update_head = update_head

        self.w[:] = decay * self.w + update_scale * xf_contiguous




    def reduction_db(self, cancelled, baseline):
        """
        Negative = reduction.
        Example: -5 dB means cancelled signal has 5 dB less power.
        """
        return 10 * np.log10(
            (np.mean(cancelled**2) + 1e-20) /
            (np.mean(baseline**2) + 1e-20)
        )


    def diagnostics(self):
        return {
            "xnorm": float(self.xnorm),
            "mavg": float(self.mavg),
            "max_step": float(self.max_step),
            "min_xnorm": float(self.min_xnorm),
            "max_control": float(self.max_control),
            "weight_norm": float(np.linalg.norm(self.w)),
        }




    def process_array(
        self,
        ref_nc,
        error_nc,
        adapt=True,
        clean_feedback=False,
        save_weight_every_s=None,
        plot=False,

        system_lag=None,
        panel_ir=None,

        poison_ref = True,

        mu=None,
        eps=None,
        leak=None,
        cancel_gain=None,
        update_sign=None,
        lag=None,
        ref_threshold=None,
        mavg_tau_ms=None,
        nonlinear=None,
    ):
        
        self.set(
            mu=mu,
            eps=eps,
            leak=leak,
            cancel_gain=cancel_gain,
            update_sign=update_sign,
            lag=lag,
            ref_threshold=ref_threshold,
            mavg_tau_ms=mavg_tau_ms,
            nonlinear=nonlinear,
        )

        ref_nc = np.asarray(ref_nc, dtype=self.dtype)
        error_nc = np.asarray(error_nc, dtype=self.dtype)

        panel_ir = self.path_ir if panel_ir is None else panel_ir

        self.adapt = bool(adapt)

        N = len(ref_nc)

        control = np.zeros(N, dtype=self.dtype)
        simulated_error = np.zeros(N, dtype=self.dtype)
        panel_output = np.zeros(N, dtype=self.dtype)

        weight_history = []

        if save_weight_every_s is not None:
            save_weight = int(save_weight_every_s * self.fs)
        else:
            save_weight = N + 1

        delayed_control_state = np.zeros(2 * self.IR_LENGTH, dtype=self.dtype)
        delayed_control_head = 0


        system_delay = int(self.lag) if system_lag is None else system_lag
        assert system_delay > 0, "System delay must be greater than 0"

        control_buf = np.zeros(system_delay, dtype=self.dtype)
        control_head = 0

       

        for n in range(N):

            delayed_control = control_buf[control_head]

            _, delayed_control_head = self._ring_push(delayed_control, delayed_control_state, delayed_control_head, self.IR_LENGTH)

            y_panel = self._ring_dot(panel_ir, delayed_control_state, delayed_control_head, self.IR_LENGTH)
            ref_feedback = self._ring_dot(self.feedback_ir, delayed_control_state, delayed_control_head, self.IR_LENGTH)


            e_n = error_nc[n] + y_panel

            poisoned_ref = ref_nc[n]
            if poison_ref:
                poisoned_ref += ref_feedback


            control_n, cleaned_ref_n = self.compute_control(poisoned_ref, clean_feedback=clean_feedback)


            control[n] = control_n
            panel_output[n] = y_panel
            simulated_error[n] = e_n

            self.process(cleaned_ref_n, e_n)


            control_buf[control_head] = control_n
            control_head += 1
            if control_head >= system_delay:
                control_head = 0


            if n % save_weight == 0:
                weight_history.append(self.w.copy())

        if save_weight_every_s is not None:
            weight_history = np.asarray(weight_history)

        db_reduction = self.reduction_db(simulated_error, error_nc)

        print(f"Simulated ANC / No ANC: {db_reduction:.2f} dB")

        if plot:
            self.plot_error_mic(error_nc, simulated_error)
            self.plot_loss_curve(error_nc, simulated_error)
            

        return (
            control,
            simulated_error,
            panel_output,
            weight_history,
        )
    






    def plot_loss_curve(
        self,
        error_nc,
        simulated_error,
        window_sec=1.0,
        hop_sec=None,
        online=None,
        save_path=None,
    ):
        """Plot error reduction over time. Positive dB means less error."""
        if hop_sec is None:
            hop_sec = window_sec

        win = int(window_sec * self.fs)
        hop = int(hop_sec * self.fs)
        if win < 1 or hop < 1:
            raise ValueError("Window and hop must each span at least one sample.")

        def reduction_over_time(cancelled):
            n = min(len(cancelled), len(error_nc))
            starts = range(0, n - win + 1, hop)
            times = (np.asarray(starts) + win / 2) / self.fs
            reductions = [
                -self.reduction_db(
                    cancelled[start:start + win],
                    error_nc[start:start + win],
                )
                for start in starts
            ]
            return times, np.asarray(reductions)

        fig, ax = plt.subplots(figsize=(6.5, 4.5))

        if online is not None:
            times, reduction = reduction_over_time(online)
            ax.plot(times, reduction, linewidth=2, label="Physical Space", color="C0")

        times, reduction = reduction_over_time(simulated_error)
        ax.plot(times, reduction, linewidth=2, label="Simulation", color="C1")

        ax.set_xlabel("Time (s)", fontsize=12)
        ax.set_ylabel("Error Reduction (dB)", fontsize=12)
        ax.set_title("Cancellation Over Time", fontsize=14, fontweight="bold")
        ax.tick_params(labelsize=10)
        ax.legend(fontsize=11, frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
        if save_path is not None:
            fig.savefig(save_path, bbox_inches="tight")

        plt.show()


    # def plot_loss_curve(self, error_nc, simulated_error, window_sec=5.0):
    #     N = min(len(simulated_error), len(error_nc))
        
    #     win = int(window_sec * self.fs)

    #     t = []
    #     loss = []
    #     for start in range(0, N - win + 1, win):
    #         end = start + win

    #         loss.append(self.reduction_db(simulated_error[start:end], error_nc[start:end]))

    #         # center of window
    #         t.append((start + win / 2) / self.fs)

    #     plt.plot(t, loss, linewidth=2, label="dB Reduction Curve", marker='o')

    #     plt.xlabel("Time (s)")
    #     plt.ylabel("Error Reduction (dB)")
    #     plt.title("Cancellation Over Time")

    #     plt.legend()
    #     plt.tight_layout()
    #     plt.show()




    def plot_error_mic(self, error_nc, simulated_error):
        plt.plot(error_nc, label="No Cancel")
        plt.plot(simulated_error, label="Simulated Cancel")

        plt.title(f"Simululated Error Mic Signal")
        plt.xlabel("Samples")
        plt.ylabel("Amplitude")
        plt.legend(loc='upper right')
        plt.grid(True)
        plt.show()
