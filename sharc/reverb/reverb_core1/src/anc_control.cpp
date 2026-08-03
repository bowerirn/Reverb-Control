#include "anc_control.h"

volatile float anc_mu            = 1e-4f;
volatile float anc_eps           = 1e-6f;
volatile float anc_leak          = 3e-7f;
volatile float anc_cancel_gain   = 0.02f;
volatile float anc_update_sign   = 1.0f;
volatile int anc_lag             = 94;
volatile bool  anc_adapt         = true;
volatile bool  anc_off         	 = false;
volatile float anc_ref_threshold = 3.0e-4f;
volatile float anc_mavg_weight   = 0.999896;
//volatile bool light = false;
