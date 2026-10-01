#include "tone.hpp"

const float TWO_PI = 6.28318530717958647692f;

// Global tone generator
Tone tone(96000.0f);

Tone::Tone(float sample_rate)
    : fs_(sample_rate),
      frequency_(500.0f),
      amplitude_(0.05f),
      phase_(0.0f),
      phase_increment_(0.0f),
      active_(false) {
    set_frequency(frequency_);
}

void Tone::set_frequency(float frequency_hz) {
    frequency_ = frequency_hz;
    phase_increment_ = TWO_PI * frequency_ / fs_;
}

void Tone::set_amplitude(float amplitude) {
    amplitude_ = amplitude;
}

float Tone::process() {
    if (!active_)
        return 0.0f;

    float output = amplitude_ * std::sin(phase_);

    phase_ += phase_increment_;

    if (phase_ >= TWO_PI)
        phase_ -= TWO_PI;

    return output;
}

void Tone::start() {
    reset();
    active_ = true;
}

void Tone::stop() {
    active_ = false;
}

bool Tone::active() const {
    return active_;
}

void Tone::reset() {
    phase_ = 0.0f;
}
