#pragma once

#include <stdint.h>

/*
 * Process one complete MIDI Control Change message.
 *
 * This may be called from the UART receive interrupt, so it must remain
 * short and nonblocking.
 */
void process_midi_control_change(
    uint8_t channel,
    uint8_t controller,
    uint8_t value
);

void anc_midi_background_loop();
