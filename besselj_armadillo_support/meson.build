armadillo_besselj_support_inc = include_directories('inc')
armadillo_dep = dependency('armadillo-code')

if not get_option('header_only')
# TODO(Antony): use "cxx.has_function" when it supports std namespace.
if not cxx.compiles(
    '''#include <cmath>

int main() {
    return int(std::cyl_bessel_j(0.0, 0.0));
}
''')
    error('Bessel function (of the first kind) does not exist in the C++ toolchain. Is this ISO C++17 compilant?')
endif

test_besselj_support_exe = executable('test_besselj_support',
    sources: [
        'tests/test-besselj.cpp',
        'tests/test-armadillo-support.cpp',
    ],
    include_directories: armadillo_besselj_support_inc,
    dependencies: [
        catch2_dep,
        armadillo_dep,
    ],
)

test('Besselj wrapper',
    test_besselj_support_exe,
    args: [
        '-r', 'tap',
    ],
    protocol: 'tap',
)
endif