from __future__ import annotations

import time
from typing import Tuple

import mido
import numpy as np
import math


class MidiProtocol:
    def __init__(self, input_name=None, output_name=None):

        input_names = mido.get_input_names()
        output_names = mido.get_output_names()

        if input_name is None:
            input_name = next(
                name for name in input_names
                if "USB2.0-MIDI" in name
            )

        if output_name is None:
            output_name = next(
                name for name in output_names
                if "USB2.0-MIDI" in name
                and "MIDIOUT2" not in name
            )

        self.input_name = input_name
        self.output_name = output_name

        print(f"\nUsing output: {self.output_name}")
        print(f"Using input:  {self.input_name}")

        try:
            self.midi_in = mido.open_input(self.input_name)
            self.midi_out = mido.open_output(self.output_name)

            self.connected = True

            while self.midi_in.poll() is not None:
                pass

        except Exception:
            print(
                "Warning: Could not open MIDI ports. "
                "Make sure the SHARC is connected and the correct ports are specified.\n"
            )

            print("MIDI outputs:")
            for i, name in enumerate(output_names):
                print(f"  {i}: {name}")

            print("\nMIDI inputs:")
            for i, name in enumerate(input_names):
                print(f"  {i}: {name}")

            self.connected = False


        


        

        # ---------------------------------------------------------------------------
        # Incoming command channels: Python -> SHARC
        # ---------------------------------------------------------------------------
        self.MIDI_ANC_MU = 0
        self.MIDI_ANC_EPS = 1
        self.MIDI_ANC_LEAK = 2
        self.MIDI_ANC_CANCEL_GAIN = 3
        self.MIDI_ANC_REF_THRESHOLD = 4
        self.MIDI_ANC_MAVG_TAU_MS = 5
        self.MIDI_ANC_UPDATE_SIGN = 6
        self.MIDI_ANC_LAG = 7
        self.MIDI_ANC_ADAPT = 8
        self.MIDI_ANC_OFF = 9
        self.MIDI_ANC_RESET = 10
        self.MIDI_ANC_GET_WEIGHTS = 11
        self.MIDI_SEED_DELTA = 12


        # ---------------------------------------------------------------------------
        # Outgoing message types: SHARC -> Python
        #
        # The SHARC puts:
        #
        #   channel bits 3-2 = message type
        #   channel bits 1-0 = top two bits of Q15 payload
        # ---------------------------------------------------------------------------
        self.MIDI_TX_WEIGHT = 0
        self.MIDI_TX_START = 1
        self.MIDI_TX_END = 2
        self.MIDI_TX_STATUS = 3


    def close(self):
        self.midi_in.close()
        self.midi_out.close()



    def _encode_u14(_, value: int) -> tuple[int, int]:
        """Split an unsigned 14-bit integer into two MIDI-safe 7-bit values."""
        if not 0 <= value <= 0x3FFF:
            raise ValueError("Value must be between 0 and 16383.")

        high = (value >> 7) & 0x7F
        low = value & 0x7F
        return high, low
    

    def _encode_small_float(_, value: float) -> tuple[int, int]:
        assert value == 0 or (1e-15 <= value < 10)

        if value == 0:
            return 0, 0

        exponent = math.ceil(-math.log10(value))
        significand = round(value * 10**exponent * 100)

        if significand == 1000:
            significand = 100
            exponent -= 1

        payload = (exponent << 10) | significand

        return (payload >> 7) & 0x7F, payload & 0x7F


    def _decode_q15_message(self, message: mido.Message) -> float:
        """
        Reconstruct the signed 16-bit Q15 payload from:

            channel bits 1-0 = payload bits 15-14
            control          = payload bits 13-7
            value            = payload bits 6-0
        """
        bits = (
            ((message.channel & 0x03) << 14)
            | ((message.control & 0x7F) << 7)
            | (message.value & 0x7F)
        )

        # Convert from two's complement
        if bits & 0x8000:
            bits -= 0x10000

        return bits / 32767.0





    def _send_control(self, channel, controller=0, value=0):
        
        """Send one Control Change message."""
        assert 0 <= channel <= 15, "MIDI channel must be between 0 and 15."
        assert 0 <= controller <= 127, "Controller must be between 0 and 127."
        assert 0 <= value <= 127, "Value must be between 0 and 127."

        self.midi_out.send(
            mido.Message(
                "control_change",
                channel=channel,
                control=controller,
                value=value,
            )
        )





    def send_small_float(self, channel, value: float):
        controller, midi_value = self._encode_small_float(value)
        self._send_control(
            channel=channel,
            controller=controller,
            value=midi_value,
        )

    def send_u14(self, channel, value: int):
        high, low = self._encode_u14(value)

        self._send_control(
            channel=channel,
            controller=high,
            value=low,
        )

    def send_bool(self, channel, value: bool):
        self._send_control(
            channel=channel,
            controller=0,
            value=int(value),
        )





    def receive_weights(self, verbose=True, timeout_s=10.0) -> tuple[np.ndarray, float]:
        expected_count = None
        weights = []

        deadline = time.monotonic() + timeout_s

        while time.monotonic() < deadline:
            message = self.midi_in.poll()

            if message is None:
                time.sleep(0.001)
                continue

            if message.type != "control_change":
                continue

            message_type = (message.channel >> 2) & 0x03

            if message_type == self.MIDI_TX_START:
                expected_count = (
                    ((message.control & 0x7F) << 7)
                    | (message.value & 0x7F)
                )

                weights.clear()
                if verbose:
                    print(f"Weight transfer started: expecting {expected_count} values")

            elif message_type == self.MIDI_TX_WEIGHT:
                weights.append(self._decode_q15_message(message))

            elif message_type == self.MIDI_TX_END:
                if len(weights) != expected_count:
                    raise RuntimeError(
                        f"Incomplete weight transfer: expected {expected_count}, "
                        f"received {len(weights)}."
                    )

                encoded_rms = self._decode_q15_message(message)
                true_rms = (encoded_rms + 1.0) / 2.0
                true_norm = true_rms * np.sqrt(expected_count)

                if verbose:
                    print(f"Weight transfer complete: received {len(weights)} values")
                return np.asarray(weights, dtype=np.float32), true_norm

        raise TimeoutError(
            f"Timed out after {timeout_s:.1f} seconds. "
            f"Received {len(weights)} weights."
        )





    def set_mu(self, mu: float):
        self.send_small_float(channel=self.MIDI_ANC_MU, value=mu)
    
    def set_leak(self, leak: float):
        self.send_small_float(channel=self.MIDI_ANC_LEAK, value=leak)

    def set_eps(self, eps: float):
        self.send_small_float(channel=self.MIDI_ANC_EPS, value=eps)

    def set_cancel_gain(self, gain: float):
        self.send_small_float(channel=self.MIDI_ANC_CANCEL_GAIN, value=gain)

    def set_ref_threshold(self, threshold: float):
        self.send_small_float(channel=self.MIDI_ANC_REF_THRESHOLD, value=threshold)   

    def set_mavg_tau_ms(self, tau_ms: int):
        self.send_u14(channel=self.MIDI_ANC_MAVG_TAU_MS, value=tau_ms)

    def set_lag(self, lag: int):
        self.send_u14(channel=self.MIDI_ANC_LAG, value=lag)
    
    def set_update_sign(self, sign: int) -> None:
        assert sign in (-1, 1), "Update sign must be -1 or +1."
        self.send_bool(channel=self.MIDI_ANC_UPDATE_SIGN, value=(sign == 1))

    def set_adapt(self, adapt: bool) -> None:
        self.send_bool(channel=self.MIDI_ANC_ADAPT, value=adapt)

    def set_off(self, off: bool) -> None:
        self.send_bool(channel=self.MIDI_ANC_OFF, value=off)



    def seed_delta(self, index: int, amplitude: float = 1.0):
        assert 0 <= index <= 2047
        assert amplitude in (-1.0, -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 1.0)

        sign = 1 if amplitude < 0 else 0

        amp_code = {
            0.25: 0,
            0.50: 1,
            0.75: 2,
            1.00: 3,
        }[abs(amplitude)]

        payload = ( (sign << 13)  |  (amp_code << 11)  |  index )

        controller = (payload >> 7) & 0x7F
        value = payload & 0x7F

        self._send_control(
            channel=self.MIDI_SEED_DELTA,
            controller=controller,
            value=value,
        )






        

    def request_reset(self) -> None:
        self._send_control(
            channel=self.MIDI_ANC_RESET,
            controller=0,
            value=1,
        )
    

    def request_weights(self, verbose=True, timeout_s=10.0) -> np.ndarray:
        self._send_control(
            channel=self.MIDI_ANC_GET_WEIGHTS,
            controller=0,
            value=1,
        )

        return self.receive_weights(verbose=verbose, timeout_s=timeout_s)

