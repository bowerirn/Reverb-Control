// #include "sine_sweep.h"
// #include <math.h>

// static const float ANC_FS = 96000.0f;
// static const float TWO_PI = 6.2831853071795864769f;

// bool anc_test_tone = false;

// static float anc_test_freq = 500.0f;
// static float anc_test_amp = 0.05f;
// static float anc_test_phase = 0.0f;


// void anc_start_test_tone(float freq, float amp)
// {
//     anc_test_freq = freq;
//     anc_test_amp = amp;
//     anc_test_phase = 0.0f;
//     anc_test_tone = true;
// }


// void anc_stop_test_tone()
// {
//     anc_test_tone = false;
//     anc_test_phase = 0.0f;
// }


// float anc_generate_test_tone()
// {
//     float y = anc_test_amp * sinf(anc_test_phase);

//     anc_test_phase += TWO_PI * anc_test_freq / ANC_FS;

//     if (anc_test_phase >= TWO_PI) {
//         anc_test_phase -= TWO_PI;
//     }

//     return y;
// }



#include "sine_sweep.h"
#include <math.h>
// #include "audio_system_config.h"


static const float FS = 96000.0f;

static const float F0 = 150.0f;
static const float F1 = 22000.0f;

static const float DURATION = 6.0f;
static const float AMPLITUDE = 0.5f;
static const float FADE_DURATION = 0.02f;

static const float TWO_PI = 6.2831853071795864769f;


SineSweep sine_sweep;


SineSweep::SineSweep()
    : running(false),
      n(0),
      N((int)(FS * DURATION)),
      fade_N((int)(FS * FADE_DURATION)),
      z(1.0f),
      q(1.0f),
      phase_scale(0.0f)
{
    float beta = logf(F1 / F0) / DURATION;

    q = expf(beta / FS);
    phase_scale = TWO_PI * F0 / beta;
}


void SineSweep::start()
{
    n = 0;
    z = 1.0f;
    running = true;
}


void SineSweep::stop()
{
    running = false;
    n = 0;
    z = 1.0f;
}


bool SineSweep::active() const
{
    return running;
}


float SineSweep::process()
{
    if (!running) {
        return 0.0f;
    }

    if (n >= N) {
        stop();
        return 0.0f;
    }

    // Same fade shape as Python make_sweep()
    float fade = 1.0f;

    if (n < fade_N) {
        fade = (float)n / (float)(fade_N - 1);
    }
    else if (n >= N - fade_N) {
        fade = (float)(N - 1 - n) / (float)(fade_N - 1);
    }

    float phase = phase_scale * (z - 1.0f);

    float y = AMPLITUDE * fade * sinf(phase);

    z *= q;
    n++;

    return y;
}
