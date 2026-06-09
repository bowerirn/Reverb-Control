import os
os.environ["SD_ENABLE_ASIO"] = "1"

import sounddevice as sd
import numpy as np

def get_asio_device(name_contains="focusrite"):
    devices = sd.query_devices()
    hostapis = sd.query_hostapis()

    for i, dev in enumerate(devices):
        hostapi_name = hostapis[dev["hostapi"]]["name"].lower()
        if "asio" in hostapi_name and name_contains.lower() in dev["name"].lower():
            return i

    raise RuntimeError(f"{name_contains} ASIO device not found")

class AudioDevice:
    def __init__(self, fs=48000, device=None):
        self.fs = fs

        self.delay_ms = 0.08983333327341825

        if device is not None:
            self.device = device
        else:
            device_id = get_asio_device()
            self.device = (device_id, device_id)
        
        sd.default.device = self.device

    @property
    def delay_samples(self):
        return self.delay_ms * self.fs

    def play(self, panel_out=None, source_out=None, dtype=np.float32, blocking=True):
        if panel_out is None:
            panel_out = np.zeros_like(source_out)
        if source_out is None:
            source_out = np.zeros_like(panel_out)

        assert len(panel_out) == len(source_out), "Input channels must have the same length"

        out = np.zeros((len(panel_out), 2), dtype=dtype)
        out[:, 0] = panel_out
        out[:, 1] = source_out

        rec = sd.playrec(
            out,
            samplerate=self.fs,
            channels=2,
            device=self.device,
            dtype=dtype,
            blocking=blocking
        )

        error_mic = rec[:, 0]
        reference_mic = rec[:, 1]

        return error_mic, reference_mic
    

    def stream_play_record(self, out0, out1, block_size=64):
        assert len(out0) == len(out1)

        N = len(out0)
        rec = np.zeros((N, 2), dtype=np.float32)
        pos = 0
        done = False

        out0 = np.asarray(out0, dtype=np.float32)
        out1 = np.asarray(out1, dtype=np.float32)

        def callback(indata, outdata, frames, time, status):
            self.delay_ms = time.outputBufferDacTime - time.inputBufferAdcTime

            nonlocal pos, done

            if status:
                print(status)

            end = pos + frames
            ncopy = max(0, min(end, N) - pos)

            outdata[:] = 0.0

            if ncopy > 0:
                outdata[:ncopy, 0] = out0[pos:pos+ncopy]
                outdata[:ncopy, 1] = out1[pos:pos+ncopy]
                rec[pos:pos+ncopy, :] = indata[:ncopy, :]

            pos = end

            if pos >= N:
                done = True

        with sd.Stream(
            samplerate=self.fs,
            blocksize=block_size,
            device=self.device,
            channels=(2, 2),
            dtype="float32",
            callback=callback,
        ):
            while not done:
                sd.sleep(50)

        return rec[:, 0], rec[:, 1]