#pragma once

#include <stddef.h>
#include "anc_control.h"



template<int N>
inline float ring_buffer_FIR(
    float new_sample,
    const float* coeffs,
    float* state_ring,
    int* head
) {
    *head = (*head == 0) ? N - 1 : *head - 1;
    state_ring[*head] = new_sample;

    float y = 0.0f;
    int len1 = N - *head;

    for (int k = 0; k < len1; k++) {
        y += coeffs[k] * state_ring[*head + k];
    }

    for (int k = len1; k < N; k++) {
        y += coeffs[k] * state_ring[k - len1];
    }

    return y;
}






template<int M, int IR_LEN, bool NLMS>
class FxLMS {
    
    public:
        FxLMS(const float* secondary_ir)
            : ir(secondary_ir),
              x_head(0),
              ir_head(0),
              xnorm(0.0f)
        {   
            reset();
        }

        void reset() {
            x_head = 0;
            ir_head = 0;
            xnorm = 0.0f;

            for (int i = 0; i < M; i++) {
                w[i] = 0.0f;
                x[i] = 0.0f;
                xf[i] = 0.0f;
            }

            for (int i = 0; i < IR_LEN; i++) {
                z[i] = 0.0f;
            }
        }

        void seed_delta(int delay, float amp) {
            if (delay >= 0 && delay < M) {
                w[delay] = amp;
            }
        }

        float process(float ref, float error_mic) {
            /*
                ref:       reference mic / clean source sample
                error_mic: current error mic sample

                returns: control sample to send to panel amp
            */

            // y = w.T * x
            float control = -anc_cancel_gain * ring_buffer_FIR<M>(ref, w, x, &x_head);

            if (!anc_adapt || anc_mu == 0.0) return control;

    

            // filtered x
            float xf_sample = ring_buffer_FIR<IR_LEN>(ref, ir, z, &ir_head);

            float step = anc_mu;
            if (NLMS) {
                float old = xf[x_head];
                xnorm += xf_sample * xf_sample - old * old;
                step = anc_mu / (anc_eps + xnorm);
            }

            float update_scale = anc_update_sign * step * error_mic;
            float decay = 1.0f - anc_leak;
            
    
            xf[x_head] = xf_sample;
            
            int update_head = x_head + anc_lag;

            while (update_head >= M) update_head -= M;
            int len1 = M - update_head;

            // w[:len1] * x[head:]
            for (int k = 0; k < len1; k++) {
                w[k] = decay * w[k] + update_scale * xf[update_head + k];
            }

            // w[len1:] * x[:head]
            for (int k = len1; k < M; k++) {
                w[k] = decay * w[k] + update_scale * xf[k - len1];
            }

            
            return control;
        }



        const float* weights() const {
            return w;
        }



    private:
        float w[M];
        float x[M];
        float xf[M];

        const float* ir;
        float z[IR_LEN];

        int x_head;
        int ir_head;

        float xnorm;

};


