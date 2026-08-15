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


template<int N>
inline float duplicated_ring_dot(
    const float* coeffs,
    const float* state,
    int head
) {
    const float* contiguous_state = state + head;

    float y = 0.0f;

    for (int k = 0; k < N; k++) {
        y += coeffs[k] * contiguous_state[k];
    }

    return y;
}


template<int N>
inline void duplicated_ring_push(
    float new_sample,
    float* state,
    int* head
) {
    int h = *head - 1;

    if (h < 0) {
        h = N - 1;
    }

    *head = h;

    state[h] = new_sample;
    state[h + N] = new_sample;
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
        float max_step;
        float min_xnorm;
        float max_control;

        FxLMS(const float* error_ir, const float* ref_ir)
            //public
            : mu(1e-4f),
              eps(1.0e-6f),
              leak(3.0e-7f),
              cancel_gain(0.02f),
              update_sign(1.0f),
              lag(78),
              adapt(false),
              ref_threshold(3.0e-4f),
              mavg_weight(0.999896f),
              max_control(0.0f),
              min_xnorm(1e30f),
              max_step(0.0f),

            //private
              path_ir(error_ir),
              feedback_ir(ref_ir),
              x_head(0),
              path_ir_head(0),
              feedback_ir_head(0),
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
            path_ir_head = 0;
            feedback_ir_head = 0;
            xnorm = 0.0f;
            mavg = 0.0f;
            update = false;

            
            max_step = 0.0f;
            min_xnorm = 1e30f;
            max_control = 0.0f;

            for (int i = 0; i < M; i++) {
                w[i] = 0.0f;
            }

            for (int i = 0; i < 2 * M; i++) {
                x[i] = 0.0f;
                xf[i] = 0.0f;
            }

            for (int i = 0; i < 2 * IR_LENGTH; i++) {
                z_path[i] = 0.0f;
            }

            for (int i = 0; i < 2 * IR_LENGTH; i++) {
                z_feedback[i] = 0.0f;
            }
        }

        void seed_delta(int delay, float amp) {
            if (delay >= 0 && delay < M) {
                w[delay] = amp;
            }
        }



        float process(float ref, float error_mic) {
            float predicted_feedback = duplicated_ring_dot<IR_LENGTH>(feedback_ir, z_feedback, feedback_ir_head);
            float cleaned_ref = ref - predicted_feedback;
            
            float control = -cancel_gain * duplicated_ring_FIR<M>(cleaned_ref, w, x, &x_head);

            duplicated_ring_push<IR_LENGTH>(control, z_feedback, &feedback_ir_head);


            if (!adapt || mu == 0.0f) {
                busy = false;
                return control;
            }


            float xf_sample = duplicated_ring_FIR<IR_LENGTH>(cleaned_ref, path_ir, z_path, &path_ir_head);

            float old = xf[x_head];
            xnorm += xf_sample * xf_sample- old * old;

            if (xnorm < 0.0f) {
                xnorm = 0.0f;
            }
            
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


            if (xnorm < min_xnorm) {
                min_xnorm = xnorm;
            }

            if (step > max_step) {
                max_step = step;
            }

            float abs_control = fabsf(control);

            if (abs_control > max_control) {
                max_control = abs_control;
            }


            busy = true;
            for (int k = 0; k < M; k++) {
                w[k] = decay * w[k] + update_scale * xf_contiguous[k];
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

        const float* path_ir;
        const float* feedback_ir;
        float z_path[2 * IR_LENGTH];
        float z_feedback[2 * IR_LENGTH];

        float block_decay;

        int x_head;
        int path_ir_head;
        int feedback_ir_head;


        float xnorm;
        float mavg;

};


