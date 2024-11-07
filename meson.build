project('MicroscPSF', 'cpp',
    license: 'MIT',
    subproject_dir: '3rdparty',
    default_options: [
        'cpp_std=c++17',
        'buildtype=debugoptimized',
    ],
)

catch2_dep = dependency('catch2-with-main')
cxx = meson.get_compiler('cpp')

subdir('besselj_armadillo_support')