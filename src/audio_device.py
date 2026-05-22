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

        if device is not None:
            self.device = device
        else:
            device_id = get_asio_device()
            self.device = (device_id, device_id)
        
        sd.default.device = self.device

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