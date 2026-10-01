#pragma once

#include <cmath>

class Tone
{
public:
    explicit Tone(float sample_rate = 96000.0f);

    void set_frequency(float frequency_hz);
    void set_amplitude(float amplitude);

    float process();

    void start();
    void stop();
    bool active() const;

    void reset();

private:
    float fs_;
    float frequency_;
    float amplitude_;

    float phase_;
    float phase_increment_;

    bool active_;
};

extern Tone tone;
