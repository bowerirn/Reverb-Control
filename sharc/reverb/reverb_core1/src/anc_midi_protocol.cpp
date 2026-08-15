#include "anc_midi_protocol.h"

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "anc_control.h"
#include "callback_midi_message.h"
#include "common/audio_system_config.h"

#include "audio_framework_selector.h"


/*
 * Each incoming Control Change channel identifies one ANC command.
 *
 * mido uses channels 0 through 15, matching the four channel bits in the
 * MIDI status byte.
 */
enum AncMidiChannel {
    MIDI_ANC_MU            = 0,
    MIDI_ANC_EPS           = 1,
    MIDI_ANC_LEAK          = 2,
    MIDI_ANC_CANCEL_GAIN   = 3,
    MIDI_ANC_REF_THRESHOLD = 4,
    MIDI_ANC_MAVG_TAU_MS   = 5,
    MIDI_ANC_UPDATE_SIGN   = 6,
    MIDI_ANC_LAG           = 7,
    MIDI_ANC_ADAPT         = 8,
    MIDI_ANC_OFF           = 9,
    MIDI_ANC_RESET         = 10,
    MIDI_ANC_GET_WEIGHTS   = 11,
    MIDI_SEED_DELTA        = 12
};


/*
 * Outgoing message type occupies the upper two MIDI-channel bits.
 *
 * The lower two channel bits carry the upper two bits of a Q15 weight.
 */
enum MidiTxType {
    MIDI_TX_WEIGHT = 0,
    MIDI_TX_START  = 1,
    MIDI_TX_END    = 2,
    MIDI_TX_STATUS = 3
};


/*
 * These are written by the MIDI receive callback and serviced later outside
 * interrupt context.
 */
static volatile bool anc_reset_requested = false;
static volatile bool anc_transfer_requested = false;


/*
 * Reconstruct a 14-bit unsigned integer from two MIDI-safe 7-bit values.
 */
static inline uint16_t decode_midi_u14(uint8_t high, uint8_t low) {
    return ((uint16_t)(high & 0x7Fu) << 7)  |  (uint16_t)(low & 0x7Fu);
}




/*
 * Decode one nonnegative floating-point value from:
 *
 *     controller: [exponent: 4 bits][significand bits 9-7]
 *     value:      [significand bits 6-0]
 *
 * The significand is stored as an integer in the range [100, 999],
 * representing values from 1.00 to 9.99. The decoded value is
 *
 *     (significand / 100) * 10^(-exponent)
 *
 * The all-zero payload (controller = 0, value = 0) represents exactly 0.0.
 */
float decode_small_float(uint8_t controller, uint8_t value) {
    const uint16_t payload =
        ((uint16_t)(controller & 0x7Fu) << 7)
        | (uint16_t)(value & 0x7Fu);

    // Reserve the all-zero representation for exactly zero.
    if (payload == 0u) {
        return 0.0f;
    }

    const uint8_t exponent = (uint8_t)((payload >> 10) & 0x0Fu);
    const uint16_t significand_code = payload & 0x03FFu;

    // Reject unused/invalid encodings.
    if (significand_code < 100u || significand_code > 999u) {
        return 0.0f;
    }

    float significand = (float)significand_code * 0.01f;

    for (uint8_t i = 0; i < exponent; i++) {
    	significand *= 0.1f;
    }

    return significand;
}


/*
 * Convert characteristic time in milliseconds into the exponential moving
 * average coefficient:
 *
 *     alpha = exp(-1 / (fs * tau))
 *
 * With tau expressed in milliseconds:
 *
 *     alpha = exp(-1000 / (fs * tau_ms))
 */
static float weight_from_tau_ms(uint16_t tau_ms) {
    if (tau_ms == 0u) {
        return 0.0f;
    }

    return expf(-1000.0f / ((float)AUDIO_SAMPLE_RATE * (float)tau_ms));
}


/*
 * Convert a weight in [-1, 1] to signed Q15.
 */
static int16_t encode_weight_q15(float weight) {
    if (weight > 1.0f)  weight = 1.0f;
    if (weight < -1.0f) weight = -1.0f;

    const float scaled = weight * 32767.0f;

    return (int16_t)(scaled >= 0.0f ? scaled + 0.5f : scaled - 0.5f);
}


/*
 * Pack one signed 16-bit Q15 value into:
 *
 *     channel:    [message type: 2 bits][weight bits 15-14]
 *     controller: [weight bits 13-7]
 *     value:      [weight bits 6-0]
 */
static void pack_q15_midi(
    int16_t quantized,
    uint8_t message_type,
    uint8_t* channel,
    uint8_t* controller,
    uint8_t* value
) {
    const uint16_t bits = (uint16_t)quantized;

    *channel = (uint8_t)(
        ((message_type & 0x03u) << 2)
        |
        ((bits >> 14) & 0x03u)
    );

    *controller = (uint8_t)((bits >> 7) & 0x7Fu);

    *value = (uint8_t)(bits & 0x7Fu);
}


/*
 * Send one encoded Q15 weight.
 *
 * midi_send_control_change() belongs to the low-level MIDI/UART layer and
 * should be declared in callback_midi_message.h.
 */
static bool send_weight(float weight) {
    const int16_t quantized = encode_weight_q15(weight);

    uint8_t channel = 0u;
    uint8_t controller = 0u;
    uint8_t value = 0u;

    pack_q15_midi(
        quantized,
        MIDI_TX_WEIGHT,
        &channel,
        &controller,
        &value
    );

    return midi_send_control_change(channel, controller, value);
}


/*
 * Start marker.
 *
 * The controller and value fields carry the 14-bit number of weights.
 */
static bool send_weight_start(uint16_t count) {
    const uint8_t channel = (uint8_t)(MIDI_TX_START << 2);
    const uint8_t high = (uint8_t)((count >> 7) & 0x7Fu);
    const uint8_t low = (uint8_t)(count & 0x7Fu);

    return midi_send_control_change(channel, high, low);
}


static bool send_weight_end(float rms) {
    const int16_t quantized = encode_weight_q15((2 * rms) - 1.0f); //[-1, 1]

    uint8_t channel = 0u;
    uint8_t controller = 0u;
    uint8_t value = 0u;

    pack_q15_midi(
        quantized,
        MIDI_TX_END,
        &channel,
        &controller,
        &value
    );

    return midi_send_control_change(channel, controller, value);
}


/*
 * Blocking transfer.
 *
 * midi_send_control_change() should return false whenever the UART transmit
 * queue cannot accept a complete message. This loop retries until each message
 * can be queued
 *
 * This function must not run inside the MIDI receive interrupt or audio ISR.
 */
static bool transfer_weights() {
    const float* w = anc.weights();
    while (!send_weight_start(FILTER_ORDER));

    float norm_sq = 0.0f;

    for (int i = 0; i < FILTER_ORDER; i++) {
        while (!send_weight(w[i]));
        norm_sq += w[i] * w[i];
    }

    float rms = sqrtf(norm_sq / FILTER_ORDER);

    while (!send_weight_end(rms));

    return true;
}


/*
 * Interpret one complete incoming MIDI Control Change message.
 *
 * This is called from interrupt context, so it only performs very short
 * parameter assignments or sets deferred request flags.
 */
void process_midi_control_change(uint8_t channel, uint8_t controller, uint8_t value) {
    
    switch (channel) {
        case MIDI_ANC_MU:
            anc.mu = decode_small_float(controller, value);
            break;

        case MIDI_ANC_EPS:
            anc.eps = decode_small_float(controller, value);
            break;

        case MIDI_ANC_LEAK:
            anc.leak = decode_small_float(controller, value);
            break;

        case MIDI_ANC_CANCEL_GAIN:
            anc.cancel_gain = decode_small_float(controller, value);
            break;

        case MIDI_ANC_REF_THRESHOLD:
            anc.ref_threshold = decode_small_float(controller, value);
            break;

        case MIDI_ANC_MAVG_TAU_MS: 
        {
            const uint16_t tau_ms = decode_midi_u14(controller, value);
            anc.mavg_weight = weight_from_tau_ms(tau_ms);
            break;
        }

        case MIDI_ANC_UPDATE_SIGN:
            anc.update_sign = (value == 0u) ? -1.0f : 1.0f;
            break;

        case MIDI_ANC_LAG:
            anc.lag = (int)decode_midi_u14(controller, value);
            break;

        case MIDI_ANC_ADAPT:
            anc.adapt = (value != 0u);
            break;

        case MIDI_ANC_OFF:
            anc_off = (value != 0u);
            break;

        case MIDI_ANC_RESET:
            anc_reset_requested = true;
            break;

        case MIDI_ANC_GET_WEIGHTS:
            anc_transfer_requested = true;
            break;

        case MIDI_SEED_DELTA:
            anc.seed_delta((int)decode_midi_u14(controller, value), 1.0f);
            break;

        default:
            break;
    }
}





void anc_midi_background_loop() {
    if (anc_transfer_requested) {
        anc_transfer_requested = false;
        anc_off = true;
        anc.adapt = false;
        transfer_weights();
        anc_off = false;
        anc.adapt = true;
    }

    if (anc_reset_requested) {
        anc_reset_requested = false;
        anc_off = true;
        anc.adapt = false;
        anc.reset();
        anc_off = false;
        anc.adapt = true;
    }
}
