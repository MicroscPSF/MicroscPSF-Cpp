microsc_psf_inc = include_directories('inc')
microsc_psf_lib = static_library('microsc-psf',
    sources: [
        'src/main.cpp',
    ],
    cpp_args: [
        '-D' + bessel_source,
    ],
    include_directories: [
        microsc_psf_inc,
        armadillo_besselj_support_inc,
    ],
    dependencies: [
        armadillo_dep,
        boost_math_dep,
    ],
)

microsc_psf_dep = declare_dependency(
    link_with: [
        microsc_psf_lib,
    ],
    include_directories: [
        microsc_psf_inc,
    ],
    dependencies: [
        armadillo_dep,
    ],
)