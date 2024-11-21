project('boost-math', 'cpp',
    version: '1.83.0',
)

boost_math_inc = include_directories('include')

boost_math_dep = declare_dependency(
    include_directories: boost_math_inc,
    compile_args: ['-DBOOST_MATH_STANDALONE'],
)