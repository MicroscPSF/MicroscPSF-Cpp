armadillo_besselj_support_inc = include_directories('inc')
armadillo_proj = subproject('armadillo-code')
armadillo_dep = armadillo_proj.get_variable('armadillo_dep')

# TODO(Antony): use "cxx.has_function" when it supports std namespace.
if get_option('use_boost')
    warning('Fallback to Boost.Math for the Bessel function implementation.')

    bessel_source = 'BOOST_HAS_BESSEL'
    boost_math_dep = dependency('boost', modules: ['math'], fallback: [
        'boost-math', 'boost_math_dep',
    ])
elif cxx.compiles(
    '''#include <cmath>
int main() {
    return int(std::cyl_bessel_j(0.0, 0.0));
}
''')
    bessel_source = 'STD_HAS_BESSEL'
    boost_math_dep = []
else
    error('Bessel function (of the first kind) does not exist in the C++ toolchain. Is this ISO C++17 compilant?')
endif

test_besselj_support_exe = executable('test_besselj_support',
    sources: [
        'tests/test-besselj.cpp',
        #'tests/test-armadillo-support.cpp',
    ],
    cpp_args: [
        '-D' + bessel_source,
    ],
    dependencies: [
        catch2_dep,
        boost_math_dep,
    ],
)

test_besselj_wrapper_exe = executable('test_besselj_armadillo_wrapper',
    sources: [
        'tests/test-armadillo-support.cpp',
    ],
    cpp_args: [
        '-D' + bessel_source,
    ],
    include_directories: [
        armadillo_besselj_support_inc,
    ],
    dependencies: [
        catch2_dep,
        armadillo_dep,
        boost_math_dep,
    ],
)

test('Besselj native support',
    test_besselj_support_exe,
    args: [
        '-r', 'tap',
    ],
    protocol: 'tap',
)

test('Besselj Armadillo wrapper',
    test_besselj_wrapper_exe,
    args: [
        '-r', 'tap',
    ],
    protocol: 'tap',
)