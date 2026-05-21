import numpy as np
from scipy.signal import chirp, fftconvolve
import sounddevice as sd


def make_sweep(fs=48000, duration=5.0, f0=20, f1=20000, fade=0.02):
    t = np.arange(int(fs * duration)) / fs

    sweep = chirp(t, f0=f0, f1=f1, t1=duration, method='logarithmic').astype(np.float32)

    # fade in/out to avoid clicks
    fade_len = int(fs * fade)
    window = np.ones_like(sweep)
    window[:fade_len] = np.linspace(0, 1, fade_len)
    window[-fade_len:] = np.linspace(1, 0, fade_len)

    sweep *= window
    sweep *= 0.5 / np.max(np.abs(sweep))

    return sweep


def estimate_ir(recorded, sweep, fs=48000, duration=5.0, f0=20, f1=20000):
    t = np.arange(len(sweep)) / fs
    k = np.log(f1/f0) / duration
    env = np.exp(t * k)

    inv = (sweep[::-1] / env).astype(np.float32)
    return fftconvolve(recorded, inv, mode='full')



def measure_ir(fs=48000, duration=5, f0=20, f1=20000, player_channel=0):
    if player_channel not in (0, 1):
        raise ValueError("player_channel must be 0 or 1")

    sweep = make_sweep(fs, duration, f0, f1)

    out = np.zeros((len(sweep), 2), dtype=np.float32)
    out[:, player_channel] = sweep

    # Play sweep, record 3 channels: [error_mic, reference_mic, loopback (no wires)]
    recording = sd.playrec(
        out,
        samplerate=fs,
        channels=3,
        dtype='float32',
        blocking=True
    )

    loopback = recording[:, 2]
    error_mic = recording[:, 0]
    reference_mic = recording[:, 1]

    loop_ir = estimate_ir(loopback, sweep, fs, duration, f0, f1)
    error_ir = estimate_ir(error_mic, sweep, fs, duration, f0, f1)
    reference_ir = estimate_ir(reference_mic, sweep, fs, duration, f0, f1)

    loop_peak = np.argmax(np.abs(loop_ir))
    error_peak = np.argmax(np.abs(error_ir))
    reference_peak = np.argmax(np.abs(reference_ir))

    hardware_latency = loop_peak
    acoustic_delay_0 = error_peak - loop_peak
    acoustic_delay_1 = reference_peak - loop_peak

    distance_m_0 = acoustic_delay_0 / fs * 343.0
    distance_m_1 = acoustic_delay_1 / fs * 343.0

    error_ir_aligned = error_ir[loop_peak:]
    reference_ir_aligned = reference_ir[loop_peak:]

    print(f"Estimated hardware latency: {hardware_latency} samples")
    print(f"Estimated distance to error mic: {distance_m_0 * 100:.2f} cm")
    print(f"Estimated distance to reference mic: {distance_m_1 * 100:.2f} cm")

    return error_ir_aligned, reference_ir_aligned, loop_ir
