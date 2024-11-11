add_languages('c')

hdf5_dep = dependency('hdf5', language: 'c', required: false)

if hdf5_dep.found()
    generate_psf_compile_args = [
        '-DARMA_USE_HDF5',
    ]
else
    generate_psf_compile_args = []
endif

generate_psf_exe = executable('generate-psf',
    sources: 'generate-psf.cpp',
    cpp_args: generate_psf_compile_args,
    dependencies: [
        microsc_psf_dep,
        hdf5_dep,
    ],
)

test('Smoke test', generate_psf_exe)