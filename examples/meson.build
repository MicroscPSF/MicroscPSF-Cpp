add_languages('c')

generate_psf_exe = executable('generate-psf',
    sources: 'generate-psf.cpp',
    dependencies: [
        microsc_psf_dep,
        dependency('hdf5', language: 'c'),
    ],
)

test('Smoke test', generate_psf_exe)