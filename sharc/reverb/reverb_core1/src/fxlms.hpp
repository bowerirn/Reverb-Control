#pragma once

#include <stddef.h>
#include "anc_control.h"



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





template<int M, int IR_LEN, int BLOCK_SIZE, bool NLMS>
class FxLMS {
    
    public:
        FxLMS(const float* secondary_ir)
            : ir(secondary_ir),
              x_head(0),
              ir_head(0),
              xnorm(0.0f),
              mavg(0.0f),
              update(false)
        {   
            reset();
        }

        void reset() {
            x_head = 0;
            ir_head = 0;
            xnorm = 0.0f;
            mavg = 0.0f;
            update = false;

            for (int i = 0; i < M; ++i) {
                w[i] = 0.0f;
            }

            for (int i = 0; i < 2 * M; ++i) {
                x[i] = 0.0f;
                xf[i] = 0.0f;
            }

            for (int i = 0; i < 2 * IR_LEN; ++i) {
                z[i] = 0.0f;
            }
        }

        void seed_delta(int delay, float amp) {
            if (delay >= 0 && delay < M) {
                w[delay] = amp;
            }
        }



        float process(float ref, float error_mic) {
        	const float mu = anc_mu;
        	const float eps = anc_eps;
        	const float leak = anc_leak;
        	const float cancel_gain = anc_cancel_gain;
        	const float update_sign = anc_update_sign;
        	const int lag = anc_lag;
        	const bool adapt = anc_adapt;
            const float mavg_weight = anc_mavg_weight;
            const float ref_threshold = anc_ref_threshold;


            float control = -cancel_gain * duplicated_ring_FIR<M>(ref, w, x, &x_head);

            if (!adapt || mu == 0.0f) return control;


            float xf_sample = duplicated_ring_FIR<IR_LEN>(ref, ir, z, &ir_head);

            float old = xf[x_head];
            xnorm += xf_sample * xf_sample- old * old;
            
            xf[x_head] = xf_sample;
            xf[x_head + M] = xf_sample;
            
            mavg = mavg_weight * mavg + (1.0f - mavg_weight) * xnorm;
            
            if (mavg < ref_threshold) {
//            	light = false;
                return control;
            }
            
//            light = true;


            
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

            for (int k = 0; k < M; k++) {
                w[k] = w[k] + update_scale * xf_contiguous[k];
            }

            return control;
        }



        const float* weights() {
            return w;
        }



    private:
        float w[M];
        float x[2 * M];
        float xf[2 * M];

        bool update;

        const float* ir;
        float z[2 * IR_LEN];

        float block_decay;

        int x_head;
        int ir_head;

        float xnorm;
        float mavg;

};


