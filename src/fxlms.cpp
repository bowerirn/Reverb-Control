// bindings.cpp
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <stdexcept>

namespace py = pybind11;

void fxlms_update(
    py::array_t<float, py::array::c_style> target,
    py::array_t<float, py::array::c_style> filtered_ref,
    py::array_t<float, py::array::c_style> x,
    py::array_t<float, py::array::c_style> y,
    py::array_t<float, py::array::c_style> e,
    py::array_t<float, py::array::c_style> w,
    float step_size,
    float eps = 1e-8f
) {
    auto target_buf = target.request();
    auto filtered_ref_buf = filtered_ref.request();
    auto x_buf = x.request();
    auto y_buf = y.request();
    auto e_buf = e.request();
    auto w_buf = w.request();

    if (target_buf.ndim != 1 || filtered_ref_buf.ndim != 1 ||
        x_buf.ndim != 1 || y_buf.ndim != 1 ||
        e_buf.ndim != 1 || w_buf.ndim != 1) {
        throw std::runtime_error("All inputs must be 1D arrays");
    }

    const int T = target_buf.shape[0];
    const int filter_order = x_buf.shape[0];

    if (filtered_ref_buf.shape[0] != T ||
        y_buf.shape[0] != T ||
        e_buf.shape[0] != T) {
        throw std::runtime_error("target, filtered_ref, y, and e must have same length");
    }

    if (w_buf.shape[0] != filter_order) {
        throw std::runtime_error("x and w must have same length");
    }

    const float* target_ptr = static_cast<float*>(target_buf.ptr);
    const float* filtered_ref_ptr = static_cast<float*>(filtered_ref_buf.ptr);

    float* x_ptr = static_cast<float*>(x_buf.ptr);
    float* y_ptr = static_cast<float*>(y_buf.ptr);
    float* e_ptr = static_cast<float*>(e_buf.ptr);
    float* w_ptr = static_cast<float*>(w_buf.ptr);


    int head = 0;

    for (int n = 0; n < T; n++) {
        head = head == 0 ? filter_order - 1 : head - 1;

        x_ptr[head] = filtered_ref_ptr[n];

        float yn = 0.0f;
        int len1 = filter_order - head;

        // dot product: x[head : end]
        for (int k = 0; k < len1; k++) {
            yn += w_ptr[k] * x_ptr[head + k];
        }

        // dot product: x[0 : head]
        for (int k = len1; k < filter_order; k++) {
            yn += w_ptr[k] * x_ptr[k - len1];
        }

        y_ptr[n] = yn;
        float err = target_ptr[n] - yn;
        e_ptr[n] = err;

        // weight update: x[head : end]
        for (int k = 0; k < len1; k++) {
            w_ptr[k] += step_size * err * x_ptr[head + k];
        }

        // weight update: x[0 : head]
        for (int k = len1; k < filter_order; k++) {
            w_ptr[k] += step_size * err * x_ptr[k - len1];
        }
    }
}



void fxnlms_update(
    py::array_t<float, py::array::c_style> target,
    py::array_t<float, py::array::c_style> filtered_ref,
    py::array_t<float, py::array::c_style> x,
    py::array_t<float, py::array::c_style> y,
    py::array_t<float, py::array::c_style> e,
    py::array_t<float, py::array::c_style> w,
    float step_size,
    float eps = 1e-8f
) {
    auto target_buf = target.request();
    auto filtered_ref_buf = filtered_ref.request();
    auto x_buf = x.request();
    auto y_buf = y.request();
    auto e_buf = e.request();
    auto w_buf = w.request();

    if (target_buf.ndim != 1 || filtered_ref_buf.ndim != 1 ||
        x_buf.ndim != 1 || y_buf.ndim != 1 ||
        e_buf.ndim != 1 || w_buf.ndim != 1) {
        throw std::runtime_error("All inputs must be 1D arrays");
    }

    const int T = target_buf.shape[0];
    const int filter_order = x_buf.shape[0];

    if (filtered_ref_buf.shape[0] != T ||
        y_buf.shape[0] != T ||
        e_buf.shape[0] != T) {
        throw std::runtime_error("target, filtered_ref, y, and e must have same length");
    }

    if (w_buf.shape[0] != filter_order) {
        throw std::runtime_error("x and w must have same length");
    }

    const float* target_ptr = static_cast<float*>(target_buf.ptr);
    const float* filtered_ref_ptr = static_cast<float*>(filtered_ref_buf.ptr);

    float* x_ptr = static_cast<float*>(x_buf.ptr);
    float* y_ptr = static_cast<float*>(y_buf.ptr);
    float* e_ptr = static_cast<float*>(e_buf.ptr);
    float* w_ptr = static_cast<float*>(w_buf.ptr);

    float xnorm = 0.0f;
    for (int k = 0; k < filter_order; k++) {
        xnorm += x_ptr[k] * x_ptr[k];
    }

    int head = 0;
    for (int n = 0; n < T; n++) {
        head = head == 0 ? filter_order - 1 : head - 1;
        float old = x_ptr[head];
        x_ptr[head] = filtered_ref_ptr[n];

        float yn = 0.0f;
        int len1 = filter_order - head;

        // dot product: x[head : end]
        for (int k = 0; k < len1; k++) {
            yn += w_ptr[k] * x_ptr[head + k];
        }

        // dot product: x[0 : head]
        for (int k = len1; k < filter_order; k++) {
            yn += w_ptr[k] * x_ptr[k - len1];
        }

        y_ptr[n] = yn;
        float err = target_ptr[n] - yn;
        e_ptr[n] = err;

        xnorm += x_ptr[head] * x_ptr[head] - old * old;
        float mu = step_size / (eps + xnorm);

        // weight update: x[head : end]
        for (int k = 0; k < len1; k++) {
            w_ptr[k] += mu * err * x_ptr[head + k];
        }

        // weight update: x[0 : head]
        for (int k = len1; k < filter_order; k++) {
            w_ptr[k] += mu * err * x_ptr[k - len1];
        }
    }
}

PYBIND11_MODULE(fxlms_ext, m) {
    m.def(
        "fxlms_update",
        &fxlms_update,
        py::arg("target"),
        py::arg("filtered_ref"),
        py::arg("x"),
        py::arg("y"),
        py::arg("e"),
        py::arg("w"),
        py::arg("step_size"),
        py::arg("eps") = 1e-8f
    );
    m.def(
        "fxnlms_update",
        &fxnlms_update,
        py::arg("target"),
        py::arg("filtered_ref"),
        py::arg("x"),
        py::arg("y"),
        py::arg("e"),
        py::arg("w"),
        py::arg("step_size"),
        py::arg("eps") = 1e-8f
    );
}
