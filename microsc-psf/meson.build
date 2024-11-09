microsc_psf_lib = static_library('microsc-psf',
    sources: [
        'src/main.cpp',
    ],
    include_directories: [
        'inc',
        armadillo_besselj_support_inc,
    ],
    dependencies: [
        armadillo_dep,
    ],
)

microsc_psf_dep = declare_dependency(
    link_with: [
        microsc_psf_lib,
    ],
    include_directories: [
        'inc',
        armadillo_besselj_support_inc,
    ],
    dependencies: [
        armadillo_dep,
    ],
)