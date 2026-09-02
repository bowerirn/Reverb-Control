#pragma once


struct Biquad {
    float b0;
    float b1;
    float b2;

    float a1;
    float a2;

    float z1;
    float z2;


    Biquad()
        : b0(0.0f),
          b1(0.0f),
          b2(0.0f),
          a1(0.0f),
          a2(0.0f),
          z1(0.0f),
          z2(0.0f)
    {
    }


    void set_coeffs(
        float b0_in,
        float b1_in,
        float b2_in,
        float a1_in,
        float a2_in
    ) {
        b0 = b0_in;
        b1 = b1_in;
        b2 = b2_in;

        a1 = a1_in;
        a2 = a2_in;
    }


    inline float process(float x) {
        float y = b0 * x + z1;

        z1 = b1 * x - a1 * y + z2;
        z2 = b2 * x - a2 * y;

        return y;
    }


    void reset() {
        z1 = 0.0f;
        z2 = 0.0f;
    }
};



class Highpass200 {

public:

    Highpass200() {

        sections[0].set_coeffs(
             0.97502896f,
            -1.95005792f,
             0.97502896f,
            -1.97485937f,
             0.97502857f
        );

        sections[1].set_coeffs(
             1.0f,
            -2.0f,
             1.0f,
            -1.98148851f,
             0.98165828f
        );

        sections[2].set_coeffs(
             1.0f,
            -2.0f,
             1.0f,
            -1.99307644f,
             0.99324720f
        );

        reset();
    }


    inline float process(float x) {

        x = sections[0].process(x);
        x = sections[1].process(x);
        x = sections[2].process(x);

        return x;
    }


    void reset() {
        for (int i = 0; i < 3; i++) {
            sections[i].reset();
        }
    }


private:

    Biquad sections[3];
};
