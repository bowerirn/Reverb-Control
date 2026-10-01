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

    float beta;
    float phase_scale;
};

extern SineSweep sine_sweep;

#endif
