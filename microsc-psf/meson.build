microsc_psf_inc = include_directories('inc')

if get_option('use_eigen')
    warning('Solving linear least squares with Eigen3, not Lapack')

    # Eigen3 is a direct competitor of Armadillo library. This may seem
    # contrived to require two nearly identicial libraries to compute PSF, but
    # Armadillo has a Matlab-like syntax that simplifies the code review process
    # a lot.
    eigen_solver_lib = static_library('eigen_solver',
        sources: 'src/linsolver.cpp',
        include_directories: include_directories('inc'),
        dependencies: [
            dependency('eigen', fallback: ['eigen', 'eigen_dep']),
            armadillo_dep,
        ],
    )

    eigen_solver_dep = declare_dependency(
        link_with: eigen_solver_lib,
        include_directories: include_directories('inc'),
        compile_args: '-DLSTSQ_USE_EIGEN',
    )
else
    eigen_solver_dep = []
endif

lapack_dep = armadillo_proj.get_variable('lapack_dep')
if lapack_dep.found()
    use_lapack_args = '-DLSTSQ_USE_LAPACK'
else
    use_lapack_args = []
endif

microsc_psf_lib = static_library('microsc-psf',
    sources: [
        'src/main.cpp',
    ],
    cpp_args: [
        '-D' + bessel_source,
        use_lapack_args,
    ],
    include_directories: [
        microsc_psf_inc,
        armadillo_besselj_support_inc,
    ],
    dependencies: [
        armadillo_dep,
        boost_math_dep,
        eigen_solver_dep,
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