python = find_program('python3')

test('exception_tests',
     python,
     args : [ files('test_exceptions.py'), dynare_preprocessor ],
     workdir : meson.current_build_dir())

if build_library.allowed()
  test('python_tests',
       python,
       args : [ '-m', 'pytest', meson.current_source_dir() ],
       env : [ 'PYTHONPATH=' + meson.project_build_root() / 'src' ],
       workdir : meson.current_source_dir(),
       timeout : 60)
endif
