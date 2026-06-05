#include <cmath>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <stdexcept>

namespace py = pybind11;

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
    float eps,
    float max_norm
) {
    auto ref_buf = ref.request();
    auto xf_in_buf = filtered_ref.request();
    auto err_buf = error_mic.request();
    auto xr_buf = x_raw.request();
    auto xf_buf = x_filt.request();
    auto cancel_buf = cancel.request();
    auto w_buf = w.request();
    auto heads_buf = heads.request();

    const int T = ref_buf.shape[0];
    const int M = w_buf.shape[0];

    const float* ref_ptr = static_cast<float*>(ref_buf.ptr);
    const float* xf_in_ptr = static_cast<float*>(xf_in_buf.ptr);
    const float* err_ptr = static_cast<float*>(err_buf.ptr);

    float* xr_ptr = static_cast<float*>(xr_buf.ptr);
    float* xf_ptr = static_cast<float*>(xf_buf.ptr);
    float* cancel_ptr = static_cast<float*>(cancel_buf.ptr);
    float* w_ptr = static_cast<float*>(w_buf.ptr);
    int* heads_ptr = static_cast<int*>(heads_buf.ptr);
    int head = heads_ptr[0];


    const float decay = 1.0f - leak;

    for (int n = 0; n < T; n++) {
        head = head == 0 ? M - 1 : head - 1;

        xr_ptr[head] = ref_ptr[n];
        xf_ptr[head] = xf_in_ptr[n];

        float y = 0.0f;
        float xnorm = 0.0f;
        int len1 = M - head;

        for (int k = 0; k < len1; k++) {
            y += w_ptr[k] * xr_ptr[head + k];
            xnorm += xf_ptr[head + k] * xf_ptr[head + k];
        }
        for (int k = len1; k < M; k++) {
            y += w_ptr[k] * xr_ptr[k - len1];
            xnorm += xf_ptr[k - len1] * xf_ptr[k - len1];
        }

        cancel_ptr[n] = -cancel_gain * y;
        float mu = step_size * err_ptr[n] / (eps + xnorm);

        for (int k = 0; k < len1; k++) {
            w_ptr[k] = w_ptr[k] * decay + mu * xf_ptr[head + k];
        }
        for (int k = len1; k < M; k++) {
            w_ptr[k] = w_ptr[k] * decay + mu * xf_ptr[k - len1];
        }
    }

    float w_norm = 0.0f;
    for (int k = 0; k < M; k++) {
        w_norm += w_ptr[k] * w_ptr[k];
    }
    w_norm = std::sqrt(w_norm);

    if (w_norm > 0.0 && w_norm > max_norm) {
        float scale = max_norm / w_norm;
        for (int k = 0; k < M; k++) {
            w_ptr[k] *= scale;
        }
        w_norm = max_norm;
    }

    heads_ptr[0] = head;
    return w_norm;
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
    float eps,
    float max_norm
) {
    auto ref_buf = ref.request();
    auto xf_in_buf = filtered_ref.request();
    auto err_buf = error_mic.request();
    auto xr_buf = x_raw.request();
    auto xf_buf = x_filt.request();
    auto cancel_buf = cancel.request();
    auto w_buf = w.request();
    auto heads_buf = heads.request();

    const int T = ref_buf.shape[0];
    const int M = w_buf.shape[0];

    const float* ref_ptr = static_cast<float*>(ref_buf.ptr);
    const float* xf_in_ptr = static_cast<float*>(xf_in_buf.ptr);
    const float* err_ptr = static_cast<float*>(err_buf.ptr);

    float* xr_ptr = static_cast<float*>(xr_buf.ptr);
    float* xf_ptr = static_cast<float*>(xf_buf.ptr);
    float* cancel_ptr = static_cast<float*>(cancel_buf.ptr);
    float* w_ptr = static_cast<float*>(w_buf.ptr);
    int* heads_ptr = static_cast<int*>(heads_buf.ptr);
    int head = heads_ptr[0];


    const float decay = 1.0f - leak;

    for (int n = 0; n < T; n++) {
        head = head == 0 ? M - 1 : head - 1;

        xr_ptr[head] = ref_ptr[n];
        xf_ptr[head] = xf_in_ptr[n];

        float y = 0.0f;
        int len1 = M - head;

        for (int k = 0; k < len1; k++) {
            y += w_ptr[k] * xr_ptr[head + k];
        }
        for (int k = len1; k < M; k++) {
            y += w_ptr[k] * xr_ptr[k - len1];
        }

        cancel_ptr[n] = -cancel_gain * y;

        float mu = step_size * err_ptr[n];

        for (int k = 0; k < len1; k++) {
            w_ptr[k] = w_ptr[k] * decay + mu * xf_ptr[head + k];
        }
        for (int k = len1; k < M; k++) {
            w_ptr[k] = w_ptr[k] * decay + mu * xf_ptr[k - len1];
        }
    }

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

    heads_ptr[0] = head;
    return w_norm;
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
        py::arg("eps") = 1e-8f,
        py::arg("max_norm") = 0.0f
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
        py::arg("eps") = 1e-8f,
        py::arg("max_norm") = 0.0f
    );
}