/*
 * Copyright (c) 2018 Analog Devices, Inc.  All rights reserved.
 *
 * These are the hooks for the MIDI / Serial processing functions.
 *
 */
#include <stdint.h>

// Define your audio system parameters in this file
#include "common/audio_system_config.h"

/**
 * UART / MIDI messages can be processed either by the ARM core or by SHARC Core 1.
 * Select which option in the audio_system_config.h file.
 */
#if defined(MIDI_UART_MANAGED_BY_SHARC1_CORE) && (MIDI_UART_MANAGED_BY_SHARC1_CORE)

// Driver for UART / MIDI functionality for Audio Project Fin
#include "drivers/bm_uart_driver/bm_uart.h"

// Event logging / error handling / functionality
#include "drivers/bm_event_logging_driver/bm_event_logging.h"

#include "callback_midi_message.h"

#include "anc_midi_protocol.h"

#include "audio_framework_selector.h"




// Create an instance of our MIDI UART driver
BM_UART midi_uart_sharc1;


bool midi_send_control_change(
    uint8_t channel,
    uint8_t controller,
    uint8_t value
) {
    channel &= 0x0Fu;
    controller &= 0x7Fu;
    value &= 0x7Fu;

    const uint8_t bytes[3] = {
        (uint8_t)(0xB0u | channel),
        controller,
        value
    };

    for (int i = 0; i < 3; i++) {
        while (uart_available_for_write(&midi_uart_sharc1) < 1u);

        if (uart_write_byte(&midi_uart_sharc1, bytes[i]) != UART_SUCCESS) {
            return false;
        }
    }

    return true;
}

/**
 * @brief Sets up MIDI on the SHARC Core 1
 *
 * @return true if successful
 */
bool midi_setup_sharc1(void)
{
    // gpio_setup(GPIO_SHARC_SAM_LED10, GPIO_OUTPUT);
	// gpio_setup(GPIO_SHARC_SAM_LED11, GPIO_OUTPUT);
	// gpio_setup(GPIO_SHARC_SAM_LED12, GPIO_OUTPUT);
    // gpio_write(GPIO_SHARC_SAM_LED10, GPIO_LOW);
    // gpio_write(GPIO_SHARC_SAM_LED11, GPIO_LOW);
	// gpio_write(GPIO_SHARC_SAM_LED12, GPIO_LOW);
    
    if (uart_initialize(
        &midi_uart_sharc1,
        UART_BAUD_RATE_MIDI,
        UART_SERIAL_8N1,
        UART_AUDIOPROJ_DEVICE_MIDI
    ) != UART_SUCCESS)
    {
        return false;
    }

    // Set our user call back for received MIDI bytes
    uart_set_rx_callback(
        &midi_uart_sharc1,
        midi_rx_callback_sharc1
    );
    
    return true;

}


/**
 * @brief Callback when new MIDI bytes arrive
 */
void midi_rx_callback_sharc1(void) {
    static uint8_t running_status = 0;
    static uint8_t controller = 0;
    static uint8_t data_count = 0;

    uint8_t byte = 0;

    while (uart_available(&midi_uart_sharc1)) {
        if (uart_read_byte(&midi_uart_sharc1, &byte) != UART_SUCCESS) {
            break;
        }

        // Status byte: bit 7 is set
        if ((byte & 0x80u) != 0u) {
            /*
             * Ignore MIDI real-time messages without interrupting
             * a partially received Control Change message.
             */
            if (byte >= 0xF8u) {
                continue;
            }

            running_status = byte;
            data_count = 0;
            continue;
        }


        // Only accept Control Change messages, 0xB0–0xBF.
        if ((running_status & 0xF0u) != 0xB0u) {
            continue;
        }

        if (data_count == 0) {
            controller = byte;
            data_count = 1;
        } else {
            const uint8_t value = byte;
            const uint8_t channel = running_status & 0x0Fu;

            data_count = 0;

            process_midi_control_change(channel, controller, value);
        }
    }
}





#endif
