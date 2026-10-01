#include "sine_sweep.h"
#include <math.h>

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
      beta(0.0f),
      phase_scale(0.0f) 
{
    beta = logf(F1 / F0) / DURATION;
    phase_scale = TWO_PI * F0 / beta;
}

void SineSweep::start() {
    n = 0;
    running = true;
}

void SineSweep::stop() {
    running = false;
    n = 0;
}

bool SineSweep::active() const {
    return running;
}

float SineSweep::process() {
    if (!running) {
        return 0.0f;
    }

    if (n >= N) {
        running = false;
        return 0.0f;
    }

    float fade_gain = 1.0f;

    if (n < fade_N) {
        fade_gain = (float)n / (float)(fade_N - 1);
    } else if (n >= N - fade_N) {
        fade_gain = (float)(N - 1 - n) / (float)(fade_N - 1);
    }

    float t = (float)n / FS;

    float z = expf(beta * t);

    float phase = phase_scale * (z - 1.0f);

    float output = AMPLITUDE * fade_gain * cosf(phase);

    n++;

    return output;
}
