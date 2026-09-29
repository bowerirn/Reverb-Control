// #ifndef SINE_SWEEP_H
// #define SINE_SWEEP_H

// // True while the test tone is replacing the normal ANC control output.
// extern bool anc_test_tone;

// // Start a sine wave at the requested frequency and digital amplitude.
// void anc_start_test_tone(float freq, float amp);

// // Stop the test tone.
// void anc_stop_test_tone();

// // Generate one output sample.
// // Call exactly once per audio sample while anc_test_tone is true.
// float anc_generate_test_tone();

// #endif




#ifndef SINE_SWEEP_H
#define SINE_SWEEP_H

class SineSweep {
public:
    SineSweep();

    void start();
    void stop();

    bool active() const;
    float process();

private:
    bool running;

    int n;
    int N;
    int fade_N;

    float z;
    float q;
    float phase_scale;
};

extern SineSweep sine_sweep;

#endif
