if have_library
  test_python = py3_inst
else
  test_python = find_program('python', 'python3', required : false)
endif

if test_python.found()
  if have_cli
    test('exception_tests',
         test_python,
         args : [ files('test_exceptions.py'), dynare_preprocessor ],
         workdir : meson.current_build_dir())
  endif

  if have_library
    test('python_tests',
         test_python,
         args : [ '-m', 'pytest', meson.current_source_dir() ],
         env : [ 'PYTHONPATH=' + meson.project_build_root() / 'src' ],
         workdir : meson.current_source_dir(),
         timeout : 60)
  endif
endif
