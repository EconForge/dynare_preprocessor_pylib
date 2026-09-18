boost_dep = dependency('boost', required : false)
if not boost_dep.found()
  # Meson's builtin boost dependency finder requires finding compiled libraries
  # even for header-only Boost usage. On macOS/Clang, compiler library search dirs
  # may not match the prefix. Fall back to finding Boost headers directly.
  fs = import('fs')
  boost_found = false
  if cpp_compiler.has_header('boost/version.hpp')
    boost_dep = declare_dependency()
    boost_found = true
  else
    py_candidate = find_program('python', 'python3', required : false)
    search_prefixes = []
    if py_candidate.found()
      env_paths = run_command(py_candidate, '-c',
        'import os; print(os.environ.get("BOOST_ROOT", "") + "\\n" + os.environ.get("CONDA_PREFIX", ""))',
        check : false).stdout().strip().split('\n')
      foreach p : env_paths
        if p != ''
          search_prefixes += [p]
        endif
      endforeach
    endif
    search_prefixes += ['/opt/homebrew', '/usr/local', '/usr']
    foreach p : search_prefixes
      if fs.is_dir(p / 'include' / 'boost')
        boost_dep = declare_dependency(include_directories : include_directories(p / 'include'))
        boost_found = true
        break
      endif
    endforeach
  endif
  if not boost_found
    error('Boost headers not found. Please install Boost or set BOOST_ROOT.')
  endif
endif
if host_machine.system() == 'emscripten'
  threads_dep = dependency('', required : false)
else
  threads_dep = dependency('threads', required : false)
endif

## Flex stuff
flex_exe = find_program('flex', 'win_flex')
# The -Ca flag comes from hitting a hard-coded size limit.
# Partial explanation: https://www.owlfolio.org/possibly-useful/flex-input-scanner-rules-are-too-complicated
# There is a Debian bug report about this: https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=642040
# Serialize Flex scanner generation to prevent win_flex temporary file collisions on Windows
tokenizer_flex = custom_target(
  'Tokenizer.cc',
  input : 'macro/Tokenizer.ll',
  output : 'Tokenizer.cc',
  command : [flex_exe, '-Ca', '-o', '@OUTPUT@', '@INPUT@']
)

dynare_flex = custom_target(
  'DynareFlex.cc',
  input : 'DynareFlex.ll',
  output : 'DynareFlex.cc',
  command : [flex_exe, '-Ca', '-o', '@OUTPUT@', '@INPUT@'],
  depends : [tokenizer_flex]
)
flex_src = [tokenizer_flex, dynare_flex]
# If FlexLexer.h is not found by the compiler, copy it to the build directory
# (adding an include directory could create problems for cross-compilation, because
# flex is typically provided by the system of the build machine)
if not cpp_compiler.has_header('FlexLexer.h')
  fs = import('fs')
  flex_candidate1 = fs.parent(fs.parent(flex_exe.full_path())) / 'include' / 'FlexLexer.h'
  flex_candidate2 = fs.parent(flex_exe.full_path()) / 'include' / 'FlexLexer.h'
  flex_candidate3 = fs.parent(flex_exe.full_path()) / 'FlexLexer.h'
  flex_candidate4 = fs.parent(fs.parent(flex_exe.full_path())) / 'share' / 'winflexbison' / 'FlexLexer.h'
  if fs.exists(flex_candidate1)
    message('FlexLexer.h copied from ' + flex_candidate1)
    flexlexer_h = fs.copyfile(flex_candidate1, 'FlexLexer.h')
  elif fs.exists(flex_candidate2)
    message('FlexLexer.h copied from ' + flex_candidate2)
    flexlexer_h = fs.copyfile(flex_candidate2, 'FlexLexer.h')
  elif fs.exists(flex_candidate3)
    message('FlexLexer.h copied from ' + flex_candidate3)
    flexlexer_h = fs.copyfile(flex_candidate3, 'FlexLexer.h')
  elif fs.exists(flex_candidate4)
    message('FlexLexer.h copied from ' + flex_candidate4)
    flexlexer_h = fs.copyfile(flex_candidate4, 'FlexLexer.h')
  else
    message('FlexLexer.h not found near flex executable')
    flexlexer_h = []
  endif
else
  flexlexer_h = []
endif

## Bison stuff
bison_exe = find_program('bison', 'win_bison')
# By convention, all our Bison source files define the api.location.file variable
# so that the generated location file follows the @BASENAME@Location.hh pattern
bison_gen = generator(bison_exe, output : ['@BASENAME@.cc', '@BASENAME@.hh', '@BASENAME@Location.hh'],
                      arguments : ['-W', '-o', '@OUTPUT0@', '@INPUT@'])
bison_src = bison_gen.process('macro/Parser.yy', 'DynareBison.yy')

preprocessor_incdir = include_directories('.', 'macro')

preprocessor_core_src = [
  'ComputingTasks.cc',
  'EquationTags.cc',
  'ModelTree.cc',
  'StaticModel.cc',
  'DynamicModel.cc',
  'NumericalConstants.cc',
  'NumericalInitialization.cc',
  'Shocks.cc',
  'SymbolTable.cc',
  'SymbolList.cc',
  'ParsingDriver.cc',
  'DataTree.cc',
  'ModFile.cc',
  'Configuration.cc',
  'Statement.cc',
  'ExprNode.cc',
  'VariableDependencyGraph.cc',
  'MacroExpandModFile.cc',
  'Bytecode.cc',
  'ExternalFunctionsTable.cc',
  'ModelEquationBlock.cc',
  'WarningConsolidation.cc',
  'SubModel.cc',
  'HeterogeneityTable.cc',
  'HeterogeneousModel.cc',
  'Utils.cc',
  'DatabaseTable.cc',
  'macro/Driver.cc',
  'macro/Environment.cc',
  'macro/Expressions.cc',
  'macro/Directives.cc'
]

if get_option('prefer_static')
  if cpp_compiler.get_id() == 'gcc' or cpp_compiler.get_id() == 'clang'
    preprocessor_link_args = [ '-static-libgcc', '-static-libstdc++' ]
    if host_machine.system() == 'windows'
      # See comment in rules for MEX files
      preprocessor_link_args += [ '-Wl,-Bstatic', '-lssp', '-Wl,-Bdynamic',
                                  '-Wl,-Bstatic,--whole-archive', '-lwinpthread', '-Wl,-Bdynamic,--no-whole-archive' ]
    endif
  else
    preprocessor_link_args = []
  endif
else
  preprocessor_link_args = []
endif

preprocessor_core_lib = static_library('dynare_core',
                                       preprocessor_core_src, flex_src, flexlexer_h, bison_src,
                                       include_directories : preprocessor_incdir,
                                       dependencies : [boost_dep, threads_dep],
                                       pic : true)

build_cli = get_option('build_cli')
build_cli_target = false
if build_cli.enabled()
  build_cli_target = true
elif build_cli.auto()
  # In auto mode, build CLI on standard targets. If cross-compiling for emscripten
  # with build_library enabled (e.g. emscripten-forge WASM side module), skip the CLI.
  if host_machine.system() != 'emscripten' or not get_option('build_library').enabled()
    build_cli_target = true
  endif
endif

have_cli = false
if build_cli_target
  dynare_preprocessor = executable('dynare-preprocessor', 'DynareMain.cc',
                                   include_directories : preprocessor_incdir,
                                   cpp_args : [ '-I' + meson.current_build_dir() / 'libdynare_core.a.p' ],
                                   dependencies : [boost_dep, threads_dep],
                                   link_with : preprocessor_core_lib,
                                   link_args : preprocessor_link_args,
                                   install : true)
  have_cli = true
endif

have_library = false
build_library = get_option('build_library')
if not build_library.disabled()
  # In Conda / Pixi on Windows, the Python executable is named 'python.exe', not 'python3.exe'.
  # Search for the interpreter running meson, 'python', and 'python3', prioritizing the one with nanobind.
  py_mod = import('python')
  py3_inst = py_mod.find_installation(required : false, pure : false)
  if not py3_inst.found() or run_command(py3_inst, '-c', 'import nanobind', check : false).returncode() != 0
    py3_candidate = py_mod.find_installation('python', required : false, pure : false)
    if py3_candidate.found() and run_command(py3_candidate, '-c', 'import nanobind', check : false).returncode() == 0
      py3_inst = py3_candidate
    else
      py3_candidate = py_mod.find_installation('python3', required : false, pure : false)
      if py3_candidate.found() and run_command(py3_candidate, '-c', 'import nanobind', check : false).returncode() == 0
        py3_inst = py3_candidate
      endif
    endif
  endif

  if py3_inst.found()
    nanobind_detect = run_command(py3_inst, '-c',
      'import nanobind, os, sys; from pathlib import Path; ' +
      'cur_dir = sys.argv[1]; src_root = sys.argv[2]; ' +
      'inc = nanobind.include_dir(); ' +
      'robin = os.path.abspath(os.path.join(inc, "..", "ext", "robin_map", "include")); ' +
      'src_dir = os.path.abspath(os.path.join(inc, "..", "src")); ' +
      'src_file = os.path.abspath(os.path.join(inc, "..", "src", "nb_combined.cpp")); ' +
      'to_meson = lambda p: (os.path.relpath(p, cur_dir) if Path(p).resolve().is_relative_to(Path(src_root).resolve()) else p).replace(\'\\\\\', \'/\'); ' +
      'print(to_meson(inc)); print(to_meson(robin)); print(to_meson(src_dir)); print(to_meson(src_file))',
      meson.current_source_dir(), meson.project_source_root(),
      check : false)

    if nanobind_detect.returncode() == 0
      nb_paths = nanobind_detect.stdout().strip().split('\n')
      nanobind_inc = include_directories(nb_paths[0], nb_paths[1], nb_paths[2])
      fs = import('fs')
      nb_combined_cpp = fs.copyfile(nb_paths[3], 'nb_combined.cpp')

      py3_dep = py3_inst.dependency(required : false)
      pylib_deps = [boost_dep, threads_dep]
      if py3_dep.found()
        pylib_deps += [py3_dep]
      endif

      dynare_pylib = py3_inst.extension_module(
        '_dynare_preprocessor',
        'DynLib.cc',
        'PythonBindings.cc',
        nb_combined_cpp,
        include_directories : [preprocessor_incdir, nanobind_inc],
        cpp_args : [ '-I' + meson.current_build_dir() / 'libdynare_core.a.p' ],
        dependencies : pylib_deps,
        link_with : preprocessor_core_lib,
        link_args : preprocessor_link_args,
        subdir : 'dynare_preprocessor',
        install : true
      )

      # Copy package source files into build directory for uninstalled test runs
      configure_file(
        input : meson.project_source_root() / 'python' / 'dynare_preprocessor' / '__init__.py',
        output : '__init__.py',
        build_subdir : 'dynare_preprocessor',
        copy : true
      )
      configure_file(
        input : meson.project_source_root() / 'python' / 'dynare_preprocessor' / 'py.typed',
        output : 'py.typed',
        build_subdir : 'dynare_preprocessor',
        copy : true
      )

      have_library = true
    elif build_library.enabled()
      error('nanobind could not be found via python3')
    endif
  endif
endif
