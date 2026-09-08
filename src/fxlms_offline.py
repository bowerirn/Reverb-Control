import numpy as np


class OfflineFxLMS:
    def __init__(
        self,
        M,
        error_ir,
        ref_ir,
        fs=96_000,
        nlms=True,
        dtype=np.float32,
    ):
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
        self.mavg_weight = dtype(0.999896)

        self.reset()

    def reset(self):
        M = self.M
        L = self.IR_LENGTH
        dtype = self.dtype

        self.x_head = 0
        self.path_ir_head = 0
        self.feedback_ir_head = 0

        self.xnorm = dtype(0.0)
        self.mavg = dtype(0.0)

        self.update = False

        # Diagnostics
        self.max_step = dtype(0.0)
        self.min_xnorm = dtype(1e30)
        self.max_control = dtype(0.0)

        # Filter weights
        self.w = np.zeros(M, dtype=dtype)

        # Duplicated ring buffers
        self.x = np.zeros(2 * M, dtype=dtype)
        self.xf = np.zeros(2 * M, dtype=dtype)

        self.z_path = np.zeros(2 * L, dtype=dtype)
        self.z_feedback = np.zeros(2 * L, dtype=dtype)

    # ------------------------------------------------------------------
    # Duplicated-ring helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _ring_push(new_sample, state, head, N):
        head -= 1

        if head < 0:
            head = N - 1

        state[head] = new_sample
        state[head + N] = new_sample

        return head

    @staticmethod
    def _ring_dot(coeffs, state, head, N):
        contiguous_state = state[head:head + N]
        return np.dot(coeffs, contiguous_state)

    @classmethod
    def _ring_fir(cls, new_sample, coeffs, state, head, N):

        head = cls._ring_push( new_sample, state, head, N)
        y = cls._ring_dot(coeffs, state, head, N)

        return y, head




    def seed_filter(self, coeffs, negate=False):
        coeffs = np.asarray(coeffs, dtype=self.dtype)

        n = min(len(coeffs), self.M)

        self.w[:n] = -coeffs[:n] if negate else coeffs[:n]
        self.w[n:] = 0.0

    def seed_delta(self, delay, amp):
        if 0 <= delay < self.M:
            self.w[delay] = self.dtype(amp)




    def compute_control(self, ref, clean_feedback=False):
        cleaned_ref = ref
        if clean_feedback:
            predicted_feedback = self._ring_dot(self.feedback_ir, self.z_feedback, self.feedback_ir_head, self.IR_LENGTH)
            cleaned_ref -= predicted_feedback

        control_raw, self.x_head = self._ring_fir(cleaned_ref, self.w, self.x, self.x_head, self.M)
        control = -self.cancel_gain * control_raw

        if clean_feedback:
            self.feedback_ir_head = self._ring_push(control, self.z_feedback, self.feedback_ir_head, self.IR_LENGTH)

        abs_control = abs(control)
        
        if abs_control > self.max_control:
            self.max_control = abs_control

        return control, cleaned_ref
    
    def process(self, cleaned_ref, error_mic):

        if not self.adapt or self.mu == 0.0:
            return

        xf_sample, self.path_ir_head = self._ring_fir(cleaned_ref, self.path_ir, self.z_path, self.path_ir_head, self.IR_LENGTH)

        old = self.xf[self.x_head]
        self.xnorm = self.xnorm + xf_sample * xf_sample - old * old

        if self.xnorm < 0.0:
            self.xnorm = 0.0

        self.xf[self.x_head] = xf_sample
        self.xf[self.x_head + self.M] = xf_sample

        self.mavg = self.mavg_weight * self.mavg + (1.0 - self.mavg_weight) * self.xnorm
    
        if self.mavg < self.ref_threshold:
            return 


        step = self.mu

        if self.NLMS:
            step = self.mu / (self.eps + self.xnorm)

        self.update = not self.update

        if not self.update:
            return 

        update_scale = self.update_sign * step * error_mic
        decay = 1.0 - self.leak

        update_head = self.x_head + self.lag

        if update_head >= self.M:
            update_head -= self.M

        xf_contiguous = self.xf[update_head:update_head + self.M]


        if self.xnorm < self.min_xnorm:
            self.min_xnorm = self.xnorm

        if step > self.max_step:
            self.max_step = step


        self.w[:] = decay * self.w + update_scale * xf_contiguous







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
    ):

        ref_nc = np.asarray(ref_nc, dtype=self.dtype)
        error_nc = np.asarray(error_nc, dtype=self.dtype)

        panel_ir = self.path_ir

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

        panel_state = np.zeros(2 * self.IR_LENGTH, dtype=self.dtype)
        panel_head = 0



        for n in range(N):

            control_n, cleaned_ref_n = self.compute_control(ref_nc[n], clean_feedback=clean_feedback)
            control[n] = control_n

            y_panel, panel_head = self._ring_fir(control_n, panel_ir, panel_state, panel_head, self.IR_LENGTH)
            panel_output[n] = y_panel

            e_n = error_nc[n] + y_panel
            simulated_error[n] = e_n

            self.process(cleaned_ref_n, e_n)

            if n % save_weight == 0:
                weight_history.append(self.w.copy())

        if save_weight_every_s is not None:
            weight_history = np.asarray(weight_history)

        return (
            control,
            simulated_error,
            panel_output,
            weight_history,
        )