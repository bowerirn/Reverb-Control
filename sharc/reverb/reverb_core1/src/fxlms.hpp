#pragma once

#include <stddef.h>
#include "math.h"


template<int N>
inline float duplicated_ring_FIR(
    float new_sample,
    const float* coeffs,
    float* state,
    int* head
) {
    // Move the head backward by one sample.
    int h = *head - 1;

    if (h < 0) {
        h = N - 1;
    }

    *head = h;

    // Maintain two identical copies of the circular buffer.
    state[h] = new_sample;
    state[h + N] = new_sample;

    // The complete N-sample history is now contiguous.
    const float* contiguous_state = state + h;

    float y = 0.0f;

    for (int k = 0; k < N; k++) {
        y += coeffs[k] * contiguous_state[k];
    }

    return y;
}





template<int M, int IR_LENGTH, bool NLMS>
class FxLMS {
    
    public:
        float mu;
        float eps;
        float leak;
        float cancel_gain;
        float update_sign;
        int lag;
        bool adapt;
        float ref_threshold;
        float mavg_weight;

        FxLMS(const float* secondary_ir)
            //public
            : mu(1e-4f),
              eps(1.0e-6f),
              leak(3.0e-7f),
              cancel_gain(0.02f),
              update_sign(1.0f),
              lag(94),
              adapt(false),
              ref_threshold(3.0e-4f),
              mavg_weight(0.999896f),

            //private
              ir(secondary_ir),
              x_head(0),
              ir_head(0),
              xnorm(0.0f),
              mavg(0.0f),
              busy(false),
              update(false)
        {   
            reset();
        }

        void reset() {
            while (busy);

            x_head = 0;
            ir_head = 0;
            xnorm = 0.0f;
            mavg = 0.0f;
            update = false;

            for (int i = 0; i < M; i++) {
                w[i] = 0.0f;
            }

            for (int i = 0; i < 2 * M; i++) {
                x[i] = 0.0f;
                xf[i] = 0.0f;
            }

            for (int i = 0; i < 2 * IR_LENGTH; i++) {
                z[i] = 0.0f;
            }
        }

        void seed_delta(int delay, float amp) {
            if (delay >= 0 && delay < M) {
                w[delay] = amp;
            }
        }



        float process(float ref, float error_mic) {

            float control = -cancel_gain * duplicated_ring_FIR<M>(ref, w, x, &x_head);

            if (!adapt || mu == 0.0f) {
                busy = false;
                return control;
            }


            float xf_sample = duplicated_ring_FIR<IR_LENGTH>(ref, ir, z, &ir_head);

            float old = xf[x_head];
            xnorm += xf_sample * xf_sample- old * old;
            
            xf[x_head] = xf_sample;
            xf[x_head + M] = xf_sample;
            
            mavg = mavg_weight * mavg + (1.0f - mavg_weight) * xnorm;
            
            if (mavg < ref_threshold) {
            	busy = false;
                return control;
            }
            

            
            float step = mu;
            
            if (NLMS) {
                step = mu / (eps + xnorm);
            }

            update = !update;

            if (!update) {
                return control;
            }

            float update_scale = update_sign * step * error_mic;
            float decay = 1.0f - leak;

            int update_head = x_head + lag;

            if (update_head >= M) {
                update_head -= M;
            }

            const float* xf_contiguous = xf + update_head;

            busy = true;
            for (int k = 0; k < M; k++) {
                w[k] = w[k] + update_scale * xf_contiguous[k];
            }
            busy = false;

            return control;
        }



        const float* weights() {
            while (busy); // This doesn't actually work you'd need a mutex but that's slow lol
            return w;
        }



    private:
        float w[M];
        float x[2 * M];
        float xf[2 * M];

        bool update;

        bool busy;

        const float* ir;
        float z[2 * IR_LENGTH];

        float block_decay;

        int x_head;
        int ir_head;

        float xnorm;
        float mavg;

};


