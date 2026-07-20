#include "anc_control.h"

volatile float anc_mu          = 1e-5f;
volatile float anc_eps         = 1e-6f;
volatile float anc_leak        = 1e-5f;
volatile float anc_cancel_gain = 0.05f;
volatile float anc_update_sign = 1.0f;
volatile int anc_lag           = 0;
volatile bool  anc_adapt       = true;
volatile bool  anc_off         = false;
