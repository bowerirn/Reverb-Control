#include <cmath>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <stdexcept>

namespace py = pybind11;

static void check_1d(const py::buffer_info& b, const char* name) {
    if (b.ndim != 1) {
        throw std::runtime_error(std::string(name) + " must be 1D");
    }
}

static float maybe_clip_w(float* w_ptr, int M, float max_norm) {
    float w_norm = 0.0f;
    for (int k = 0; k < M; k++) {
        w_norm += w_ptr[k] * w_ptr[k];
    }
    w_norm = std::sqrt(w_norm);

    if (max_norm > 0.0f && w_norm > max_norm) {
        float scale = max_norm / w_norm;
        for (int k = 0; k < M; k++) {
            w_ptr[k] *= scale;
        }
        w_norm = max_norm;
    }

    return w_norm;
}

float nlms_realtime_update(
    py::array_t<float, py::array::c_style> ref,
    py::array_t<float, py::array::c_style> filtered_ref,
    py::array_t<float, py::array::c_style> error_mic,
    py::array_t<float, py::array::c_style> x_raw,
    py::array_t<float, py::array::c_style> x_filt,
    py::array_t<float, py::array::c_style> cancel,
    py::array_t<float, py::array::c_style> w,
    py::array_t<int, py::array::c_style> heads,
    float step_size,
    float leak,
    float cancel_gain,
    float w_update_sign,
    float eps,
    float max_norm,
    int adapt_delay,
    int control_delay
) {
    auto ref_buf = ref.request();
    auto xf_in_buf = filtered_ref.request();
    auto err_buf = error_mic.request();
    auto x_buf = x_raw.request();
    auto xf_buf = x_filt.request();
    auto cancel_buf = cancel.request();
    auto w_buf = w.request();
    auto heads_buf = heads.request();

    check_1d(ref_buf, "ref");
    check_1d(xf_in_buf, "filtered_ref");
    check_1d(err_buf, "error_mic");
    check_1d(x_buf, "x_raw");
    check_1d(xf_buf, "x_filt");
    check_1d(cancel_buf, "cancel");
    check_1d(w_buf, "w");
    check_1d(heads_buf, "heads");

    const int T = static_cast<int>(ref_buf.shape[0]);
    const int M = static_cast<int>(w_buf.shape[0]);
    const int Mf = static_cast<int>(xf_buf.shape[0]);
    const int Mr = static_cast<int>(x_buf.shape[0]);


    if (xf_in_buf.shape[0] != T || err_buf.shape[0] != T || cancel_buf.shape[0] != T) {
        throw std::runtime_error("ref, filtered_ref, error_mic, and cancel must have same length");
    }
    if (Mf < M + adapt_delay) {
        throw std::runtime_error("x_filt length must be at least M + adapt_delay");
    }
    if (Mr < M + control_delay) {
        throw std::runtime_error("x length must be at least M + control_delay");
    }
    if (heads_buf.shape[0] < 2) {
        throw std::runtime_error("heads must have length at least 2");
    }

    const float* ref_ptr = static_cast<float*>(ref_buf.ptr);
    const float* xf_in_ptr = static_cast<float*>(xf_in_buf.ptr);
    const float* err_ptr = static_cast<float*>(err_buf.ptr);

    float* x_ptr = static_cast<float*>(x_buf.ptr);
    float* xf_ptr = static_cast<float*>(xf_buf.ptr);
    float* cancel_ptr = static_cast<float*>(cancel_buf.ptr);
    float* w_ptr = static_cast<float*>(w_buf.ptr);
    int* heads_ptr = static_cast<int*>(heads_buf.ptr);

    int raw_head = heads_ptr[0];
    int filt_head = heads_ptr[1];

    const float decay = 1.0f - leak;

    for (int n = 0; n < T; n++) {
        raw_head = (raw_head == 0) ? Mr - 1 : raw_head - 1;
        filt_head = (filt_head == 0) ? Mf - 1 : filt_head - 1;

        x_ptr[raw_head] = ref_ptr[n];
        xf_ptr[filt_head] = xf_in_ptr[n];

        int output_head = raw_head + control_delay;
        while (output_head >= Mr) output_head -= Mr;

        float y = 0.0f;
        for (int k = 0; k < M; k++) {
            int idx = output_head + k;
            if (idx >= Mr) idx -= Mr;
            y += w_ptr[k] * x_ptr[idx];
        }

        cancel_ptr[n] = -cancel_gain * y;

        int update_head = filt_head + adapt_delay;
        while (update_head >= Mf) update_head -= Mf;

        float xnorm = 0.0f;
        for (int k = 0; k < M; k++) {
            int idx = update_head + k;
            if (idx >= Mf) idx -= Mf;
            xnorm += xf_ptr[idx] * xf_ptr[idx];
        }

        float mu = step_size * err_ptr[n] / (eps + xnorm);

        for (int k = 0; k < M; k++) {
            int idx = update_head + k;
            if (idx >= Mf) idx -= Mf;
            w_ptr[k] = w_ptr[k] * decay + w_update_sign * mu * xf_ptr[idx];
        }
    }

    heads_ptr[0] = raw_head;
    heads_ptr[1] = filt_head;

    return maybe_clip_w(w_ptr, M, max_norm);
}

float lms_realtime_update(
    py::array_t<float, py::array::c_style> ref,
    py::array_t<float, py::array::c_style> filtered_ref,
    py::array_t<float, py::array::c_style> error_mic,
    py::array_t<float, py::array::c_style> x_raw,
    py::array_t<float, py::array::c_style> x_filt,
    py::array_t<float, py::array::c_style> cancel,
    py::array_t<float, py::array::c_style> w,
    py::array_t<int, py::array::c_style> heads,
    float step_size,
    float leak,
    float cancel_gain,
    float w_update_sign,
    float eps,
    float max_norm,
    int adapt_delay,
    int control_delay
) {
    auto ref_buf = ref.request();
    auto xf_in_buf = filtered_ref.request();
    auto err_buf = error_mic.request();
    auto x_buf = x_raw.request();
    auto xf_buf = x_filt.request();
    auto cancel_buf = cancel.request();
    auto w_buf = w.request();
    auto heads_buf = heads.request();

    check_1d(ref_buf, "ref");
    check_1d(xf_in_buf, "filtered_ref");
    check_1d(err_buf, "error_mic");
    check_1d(x_buf, "x_raw");
    check_1d(xf_buf, "x_filt");
    check_1d(cancel_buf, "cancel");
    check_1d(w_buf, "w");
    check_1d(heads_buf, "heads");

    const int T = static_cast<int>(ref_buf.shape[0]);
    const int M = static_cast<int>(w_buf.shape[0]);
    const int Mf = static_cast<int>(xf_buf.shape[0]);
    const int Mr = static_cast<int>(x_buf.shape[0]);

    if (xf_in_buf.shape[0] != T || err_buf.shape[0] != T || cancel_buf.shape[0] != T) {
        throw std::runtime_error("ref, filtered_ref, error_mic, and cancel must have same length");
    }
    if (Mf < M + adapt_delay) {
        throw std::runtime_error("x_filt length must be at least M + adapt_delay");
    }
    if (Mr < M + control_delay) {
        throw std::runtime_error("x length must be at least M + control_delay");
    }
    if (heads_buf.shape[0] < 2) {
        throw std::runtime_error("heads must have length at least 2");
    }

    const float* ref_ptr = static_cast<float*>(ref_buf.ptr);
    const float* xf_in_ptr = static_cast<float*>(xf_in_buf.ptr);
    const float* err_ptr = static_cast<float*>(err_buf.ptr);

    float* x_ptr = static_cast<float*>(x_buf.ptr);
    float* xf_ptr = static_cast<float*>(xf_buf.ptr);
    float* cancel_ptr = static_cast<float*>(cancel_buf.ptr);
    float* w_ptr = static_cast<float*>(w_buf.ptr);
    int* heads_ptr = static_cast<int*>(heads_buf.ptr);

    int raw_head = heads_ptr[0];
    int filt_head = heads_ptr[1];

    const float decay = 1.0f - leak;

    for (int n = 0; n < T; n++) {
        raw_head = (raw_head == 0) ? Mr - 1 : raw_head - 1;
        filt_head = (filt_head == 0) ? Mf - 1 : filt_head - 1;

        x_ptr[raw_head] = ref_ptr[n];
        xf_ptr[filt_head] = xf_in_ptr[n];

        int output_head = raw_head + control_delay;
        while (output_head >= Mr) output_head -= Mr;

        float y = 0.0f;
        for (int k = 0; k < M; k++) {
            int idx = output_head + k;
            if (idx >= Mr) idx -= Mr;
            y += w_ptr[k] * x_ptr[idx];
        }

        cancel_ptr[n] = -cancel_gain * y;

        int update_head = filt_head + adapt_delay;
        while (update_head >= Mf) update_head -= Mf;

        float mu = step_size * err_ptr[n];

        for (int k = 0; k < M; k++) {
            int idx = update_head + k;
            if (idx >= Mf) idx -= Mf;
            w_ptr[k] = w_ptr[k] * decay + w_update_sign * mu * xf_ptr[idx];
        }
    }

    heads_ptr[0] = raw_head;
    heads_ptr[1] = filt_head;

    return maybe_clip_w(w_ptr, M, max_norm);
}

PYBIND11_MODULE(lms_rt_ext, m) {
    m.def(
        "nlms_realtime_update",
        &nlms_realtime_update,
        py::arg("ref"),
        py::arg("filtered_ref"),
        py::arg("error_mic"),
        py::arg("x_raw"),
        py::arg("x_filt"),
        py::arg("cancel"),
        py::arg("w"),
        py::arg("heads"),
        py::arg("step_size"),
        py::arg("leak") = 0.0f,
        py::arg("cancel_gain") = 1.0f,
        py::arg("w_update_sign") = 1.0f,
        py::arg("eps") = 1e-8f,
        py::arg("max_norm") = 0.0f,
        py::arg("adapt_delay") = 0,
        py::arg("control_delay") = 0
    );

    m.def(
        "lms_realtime_update",
        &lms_realtime_update,
        py::arg("ref"),
        py::arg("filtered_ref"),
        py::arg("error_mic"),
        py::arg("x_raw"),
        py::arg("x_filt"),
        py::arg("cancel"),
        py::arg("w"),
        py::arg("heads"),
        py::arg("step_size"),
        py::arg("leak") = 0.0f,
        py::arg("cancel_gain") = 1.0f,
        py::arg("w_update_sign") = 1.0f,
        py::arg("eps") = 1e-8f,
        py::arg("max_norm") = 0.0f,
        py::arg("adapt_delay") = 0,
        py::arg("control_delay") = 0
    );
}