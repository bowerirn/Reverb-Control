import numpy as np
from scipy.signal import chirp, fftconvolve
from .audio_device import AudioDevice


def make_sweep(fs=48000, duration=6.0, f0=150, f1=22000, fade=0.02):
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


def estimate_ir(recorded, sweep, fs=48000, f0=150, f1=22000):
    duration = len(sweep) / fs
    t = np.arange(len(sweep)) / fs
    k = np.log(f1 / f0) / duration

    env = np.exp(-t * k)
    inv = (sweep[::-1] * env).astype(np.float32)

    ir = fftconvolve(recorded, inv, mode="full")

    # normalize so sweep convolved with inverse has peak 1
    auto = fftconvolve(sweep, inv, mode="full")
    norm = np.max(np.abs(auto))

    return ir / norm



def measure_ir(ad: AudioDevice, duration=6.0, f0=150, f1=22000, player_channel=0, silence=4.0):
    if player_channel not in (0, 1):
        raise ValueError("player_channel must be 0 or 1")

    sweep = make_sweep(ad.fs, duration, f0, f1)

    decay_silence = np.zeros(int(ad.fs * silence), dtype=np.float32)  # 2 sec decay
    play_signal = np.concatenate([sweep, decay_silence])

    out = np.zeros((len(play_signal), 2), dtype=np.float32)
    out[:len(sweep), player_channel] = sweep

    error_mic, ref_mic = ad.play(out[:, 0], out[:, 1])

    error_ir = estimate_ir(error_mic, sweep, ad.fs, f0, f1)
    ref_ir = estimate_ir(ref_mic, sweep, ad.fs, f0, f1)


    return error_ir, ref_ir




def measure_ir_stream(
    ad: AudioDevice,
    duration=6.0,
    f0=150,
    f1=22000,
    player_channel=0,
    silence=4.0,
    block_size=64,
):
    if player_channel not in (0, 1):
        raise ValueError("player_channel must be 0 or 1")

    sweep = make_sweep(ad.fs, duration, f0, f1)

    decay_silence = np.zeros(int(ad.fs * silence), dtype=np.float32)
    play_signal = np.concatenate([sweep, decay_silence])

    out0 = np.zeros(len(play_signal), dtype=np.float32)
    out1 = np.zeros(len(play_signal), dtype=np.float32)

    if player_channel == 0:
        out0[:len(sweep)] = sweep
    else:
        out1[:len(sweep)] = sweep

    error_mic, ref_mic = ad.stream_play_record(
        out0, out1, block_size=block_size,
    )

    error_ir = estimate_ir(error_mic, sweep, ad.fs, f0, f1)
    ref_ir = estimate_ir(ref_mic, sweep, ad.fs, f0, f1)

    return error_ir, ref_ir



def align_ir_by_distance(ir, distance_cm, ir_len=128, fs=48000):
    sound_speed = 343.0

    acoustic_delay = int(round((distance_cm / 100.0) / sound_speed * fs))

    peak = np.argmax(np.abs(ir))

    start = peak - acoustic_delay
    if start < 0:
        start = 0

    end = start + ir_len

    if end > len(ir):
        out = np.zeros(ir_len, dtype=ir.dtype)
        available = ir[start:]
        out[:len(available)] = available
        return out

    return ir[start:end]
