from . import fxlms_ext
import numpy as np
from scipy.signal import lfilter
import matplotlib.pyplot as plt



class FxLMS:
    def __init__(
        self, 
        ref, 
        target, 
        ir, 
        filter_order, 
        step_fn, 
        use_norm = True, 
        eps = 1e-8
    ):
        self.ref = np.ascontiguousarray(ref.squeeze(), dtype=np.float32)
        self.target = np.ascontiguousarray(target.squeeze(), dtype=np.float32)
        self.ir = np.ascontiguousarray(ir.squeeze(), dtype=np.float32)
        self.filter_order = filter_order
        self.use_norm = use_norm
        self.eps = eps

        if isinstance(step_fn, (int, float)):
            self.step_fn = lambda i: step_fn
        else: self.step_fn = step_fn

        self.ref_filt = np.ascontiguousarray(
            lfilter(self.ir, np.array([1.0], dtype=np.float32), self.ref),
            dtype=np.float32
        )
        self.x = np.zeros(self.filter_order, dtype=np.float32)
        self.w = np.zeros(self.filter_order, dtype=np.float32)
        self.y = np.zeros(len(self.target), dtype=np.float32)
        self.e = np.zeros(len(self.target), dtype=np.float32)

        self.update = fxlms_ext.fxnlms_update if use_norm else fxlms_ext.fxlms_update
        self.mses = []

    @property
    def pred(self):
        return lfilter(self.w, [1.0], self.ref)

    def reset(self):
        self.x.fill(0)
        self.w.fill(0)
        self.y.fill(0)
        self.mses = []


    def learn(self, n_iter, reset=True, verbose=False, plot_result=True, plot_mse=True):
        if reset:
            self.reset()

        for i in range(n_iter):
            if verbose:
                print(f"Iteration {i + 1}/{n_iter}")

            self.update(
                self.target,
                self.ref_filt,
                self.x,
                self.y,
                self.e,
                self.w,
                np.float32(self.step_fn(i + 1)),
                np.float32(self.eps),
            )

            y_eval = lfilter(self.w, [1.0], self.ref_filt)
            e_eval = self.target - y_eval
            mse = np.mean(e_eval ** 2)
            self.mses.append(mse)

        print(f'Final MSE: {self.mses[-1]:.3e}')
        print(f'Lowest MSE: {min(self.mses):.3e}')

        if plot_result:
            self.results_plot()

        if plot_mse:
            self.plot_mses()


    


    def results_plot(self):
        pred = lfilter(self.w, [1.0], self.ref_filt)

        plt.figure()
        plt.plot(self.target, label="ground truth", alpha=0.7)
        plt.plot(pred, label="learned prediction", alpha=0.7)
        plt.legend()
        plt.title('FxLMS Results')
        plt.xlabel('Samples')
        plt.ylabel('Amplitude')
        plt.show()


    def plot_mses(self):
        mses = np.asarray(self.mses)

        plt.figure(figsize=(8, 4))
        plt.plot(mses)
        plt.xlabel("Iteration")
        plt.ylabel("MSE")
        plt.title("MSE vs Iteration")
        plt.grid(True)
        plt.tight_layout()
        plt.show()