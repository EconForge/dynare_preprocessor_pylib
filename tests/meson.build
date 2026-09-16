python = find_program('python3')

test('exception_tests',
     python,
     args : [ files('test_exceptions.py'), dynare_preprocessor ],
     workdir : meson.current_build_dir())
