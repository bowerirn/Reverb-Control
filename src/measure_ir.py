import os
os.environ["SD_ENABLE_ASIO"] = "1"

import numpy as np
from scipy.signal import chirp, fftconvolve, correlate
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


def estimate_ir(recorded, sweep, fs=48000, f0=20, f1=20000):
    duration = len(sweep) / fs
    t = np.arange(len(sweep)) / fs
    k = np.log(f1 / f0) / duration

    # decay envelope applied to reversed sweep
    env = np.exp(-t * k)

    inv = (sweep[::-1] * env).astype(np.float32)
    ir = fftconvolve(recorded, inv, mode="full")

    return ir



def compute_latency(latency_sweep, sweep_rec, fs=48000, f0=24, f1=2400):
    # MATLAB-like loopback quality check
    xc = correlate(sweep_rec, latency_sweep, mode="full")
    peak = np.max(np.abs(xc))
    quality = peak / (np.sqrt(np.mean(xc**2)) + 1e-12)

    if peak < 0.01 or quality < 10:
        raise RuntimeError(
            "Loopback signal is too weak/noisy. Check routing or physical loopback cable."
        )

    loop_ir = estimate_ir(sweep_rec, latency_sweep, fs, f0, f1)
    latency = int(np.argmax(np.abs(loop_ir)) - (len(latency_sweep) - 1))

    return latency, loop_ir



def measure_ir(fs=48000, duration=5, f0=20, f1=20000, latency_duration=0.5, decay_duration=2.0, player_channel=0, latency_channel=0):

    if player_channel not in (0, 1):
        raise ValueError("player_channel must be 0 or 1")
    
    latency_f0 = fs / 2000
    latency_f1 = fs / 20

    latency_sweep = 0.5 * make_sweep(
        fs=fs,
        duration=latency_duration,
        f0=latency_f0,
        f1=latency_f1,
    )


    sweep = make_sweep(fs, duration, f0, f1)

    pre_silence = np.zeros(max(int(0.1 * fs), 1024), dtype=np.float32)
    decay = np.zeros(int(decay_duration * fs), dtype=np.float32)

    # Timeline:
    # silence -> latency_sweep -> main_sweep -> decay silence
    total_len = len(pre_silence) + len(latency_sweep) + len(sweep) + len(decay)
    out = np.zeros((total_len, 2), dtype=np.float32)

    latency_start = len(pre_silence)
    main_start = latency_start + len(latency_sweep)

    out[latency_start:latency_start + len(latency_sweep), player_channel] = latency_sweep
    out[main_start:main_start + len(sweep), player_channel] = sweep


    # Play sweep, record 3 channels: [error_mic, reference_mic, loopback (no wires)]
    recording = sd.playrec(
        out,
        samplerate=fs,
        channels=3,
        dtype='float32',
        blocking=True
    )

    latency_rec = recording[
        latency_start:latency_start + len(latency_sweep),
        latency_channel,
    ]

    error_arrival_samples, error_arrival_ir = compute_latency(
        latency_sweep,
        latency_rec,
        fs=fs,
        f0=latency_f0,
        f1=latency_f1,
    )

    main_rec = recording[main_start:, :]

    error_mic = main_rec[:, 0]
    reference_mic = main_rec[:, 1]

    error_ir = estimate_ir(error_mic, sweep, fs, f0, f1)
    ref_ir = estimate_ir(reference_mic, sweep, fs, f0, f1)


    error_ir_aligned = error_ir[error_arrival_samples:]
    ref_ir_aligned = ref_ir[error_arrival_samples:]


    print(f"Error mic arrival delay: {error_arrival_samples} samples")
    print(f"Error mic arrival delay: {error_arrival_samples / fs * 1000:.3f} ms")

    return {
        "error_ir": error_ir,
        "ref_ir": ref_ir,
        "error_ir_aligned": error_ir_aligned,
        "ref_ir_aligned": ref_ir_aligned,
        "error_arrival_ir": error_arrival_ir,
        "latency_samples": error_arrival_samples,
        "recording": recording,
        "main_sweep": sweep,
        "latency_sweep": latency_sweep,
    }
