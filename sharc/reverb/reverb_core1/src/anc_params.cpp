#include "anc_params.h"

volatile float anc_mu            = 1e-4f;
volatile float anc_eps           = 1.0e-6f;
volatile float anc_leak          = 3.0e-7f;
volatile float anc_cancel_gain   = 0.02f;
volatile float anc_update_sign   = 1.0f;
volatile int   anc_lag           = 94;
volatile bool  anc_adapt         = false;
volatile bool  anc_off           = true;
volatile float anc_ref_threshold = 3.0e-4f;
volatile float anc_mavg_weight   = 0.999896f;
