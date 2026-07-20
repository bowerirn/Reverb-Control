import os
os.environ["SD_ENABLE_ASIO"] = "1"

import sounddevice as sd
import numpy as np
import gc

def get_device(name_contains="focusrite"):
    if isinstance(name_contains, int):
        return name_contains
    
    devices = sd.query_devices()
    # hostapis = sd.query_hostapis()

    for i, dev in enumerate(devices):
        # hostapi_name = hostapis[dev["hostapi"]]["name"].lower()
        if name_contains.lower() in dev["name"].lower():
            return i

    raise RuntimeError(f"{name_contains} device not found")

class AudioDevice:
    def __init__(self, fs=48_000, in_device=None, out_device=None):

        self.system_delay_ms = 89.83333327341825

        self.available_samplerates = [44_100, 48_000, 88_200, 96_000, 176_400, 192_000]
        self.available_buffer_sizes = [16, 32, 48, 64, 96, 128, 160, 192, 256, 512, 1024]
        self.device_latencies_s = [.0049, .0055, .0062, .0079, .0092, .0145, .0179, .0192, .0259, .0465, .0899]

        self.fs = fs
        if self.fs not in self.available_samplerates:
            self.fs = 48_000
            print(f"Sample rate of {fs} is not available, defaulting to 48kHz.\nAvailable sample rates are: [44.1kHz, 48kHz, 88.2kHz, 96kHz, 176.4kHz, 192KhZ]")


        self.per_buffer_latency_s = {
            buffer_size: latency_s
            for buffer_size, latency_s in zip(
                self.available_buffer_sizes,
                self.device_latencies_s
            )
        }

        
        if in_device is None:
            in_device = 'focusrite'
        if out_device is None:
            out_device = 'focusrite'
        else:
            in_device_id = get_device(in_device)
            out_device_id = get_device(out_device)
            self.device = (in_device_id, out_device_id)
        
        sd.default.device = self.device

    @property
    def delay_samples(self):
        return self.delay_ms * (self.fs / 1000)

    def play(self, left=None, right=None, dtype=np.float32, blocking=True):
        if left is None and right is None:
            return
        
        if left is None:
            left = np.zeros_like(right)
        elif right is None:
            right = np.zeros_like(left)

        assert len(left) == len(right), "Input channels must have the same length"

        out = np.zeros((len(left), 2), dtype=dtype)
        out[:, 0] = left
        out[:, 1] = right

        rec = sd.playrec(
            out,
            samplerate=self.fs,
            channels=2,
            device=self.device,
            dtype=dtype,
            blocking=blocking,
        )

        error_mic = rec[:, 0]
        reference_mic = rec[:, 1]

        return error_mic, reference_mic
    

    def stream(self, callback, block_size=192, dtype="float32"):
        assert block_size in self.available_buffer_sizes, f"Device does not support block size of {block_size}. Please use one of: {str(self.available_buffer_sizes)}"

        with sd.Stream(
            samplerate=self.fs,
            blocksize=block_size,
            device=self.device,
            channels=(2, 2),
            dtype=dtype,
            latency=(0.0, 0.0),
            callback=callback,
        ) as stream:
            while stream.active:
                sd.sleep(100)
        
        gc.collect()

    def stop_stream(self):
        raise sd.CallbackStop()