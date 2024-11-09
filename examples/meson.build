generate_psf_exe = executable('generate-psf',
    sources: 'generate-psf.cpp',
    dependencies: [
        microsc_psf_dep,
    ],
)

test('Smoke test', generate_psf_exe)