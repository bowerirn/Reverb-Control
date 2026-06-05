from setuptools import setup, Extension
import pybind11

ext_modules = [
    Extension(
        "src.fxlms_ext",
        ["src/fxlms.cpp"],
        include_dirs=[pybind11.get_include()],
        language="c++",
        extra_compile_args=["-O3", "-std=c++17", "-march=native"],
    ),
    Extension(
        "src.lms_rt_ext",
        ["src/rt_lms.cpp"],
        include_dirs=[pybind11.get_include()],
        language="c++",
        extra_compile_args=["-O3", "-std=c++17", "-march=native"],
    )
]

setup(
    ext_modules=ext_modules,
)
