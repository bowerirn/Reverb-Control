/*
 * Copyright (c) 2018 Analog Devices, Inc.  All rights reserved.
 */
#ifndef _CALLBACK_MIDI_MESSAGE_H
#define _CALLBACK_MIDI_MESSAGE_H

// Define your audio system parameters in this file
#include "common/audio_system_config.h"

// UART functionality for MIDI driver on Audio Project Fin
#include "drivers/bm_uart_driver/bm_uart.h"

#if (MIDI_UART_MANAGED_BY_SHARC1_CORE)

#ifdef __cplusplus
extern "C" {
#endif

extern BM_UART MIDI_UART_SHARC1;

// Sets up the MIDI
bool midi_setup_sharc1(void);

// MIDI callback for received event
void midi_rx_callback_sharc1(void);

bool midi_send_control_change(uint8_t channel, uint8_t controller, uint8_t value);

#ifdef __cplusplus
}
#endif

#endif

#endif     // _CALLBACK_MIDI_MESSAGE_H
