project('MicroscPSF', 'cpp',
    license: 'MIT',
    license_files: [
        'LICENSE.txt',
        '3rdparty/armadillo-code/LICENSE.txt',
    ],
    subproject_dir: '3rdparty',
    default_options: [
        'cpp_std=c++17',
        'buildtype=debugoptimized',
    ],
)

catch2_dep = dependency('catch2-with-main')
cxx = meson.get_compiler('cpp')

subdir('besselj_armadillo_support')
subdir('microsc-psf')
subdir('examples')