# TODO for reaching feature parity with current autotools-based build system:
# - Configuration option to disable documentation
# - Add -Wold-style-cast flag except when building flex-generated files

project('dynare-preprocessor', 'cpp',
        version : '6-unstable',
        default_options : [ 'cpp_std=gnu++20', 'warning_level=2' ])

add_global_arguments('-DPACKAGE_VERSION="' + meson.project_version() + '"', language : 'cpp')

boost_dep = dependency('boost')

## Flex stuff
flex_exe = find_program('flex')
# The -Ca flag comes from hitting a hard-coded size limit.
# Partial explanation: https://www.owlfolio.org/possibly-useful/flex-input-scanner-rules-are-too-complicated
# There is a Debian bug report about this: https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=642040
flex_gen = generator(flex_exe, output : '@BASENAME@.cc',
                     arguments : ['-Ca', '-o', '@OUTPUT@', '@INPUT@'])
flex_src = flex_gen.process('src/macro/Tokenizer.ll', 'src/DynareFlex.ll')
# If FlexLexer.h is not found by the compiler, copy it to the build directory
# (adding an include directory could create problems for cross-compilation, because
# flex is typically provided by the system of the build machine)
compiler = meson.get_compiler('cpp')
if not compiler.has_header('FlexLexer.h')
  message('FlexLexer.h cannot be found by the compiler, it will be manually copied to the build directory')
  fs = import('fs')
  flexlexer_h = fs.copyfile(fs.parent(fs.parent(flex_exe.full_path())) / 'include' / 'FlexLexer.h', 'FlexLexer.h')
else
  flexlexer_h = []
endif

## Bison stuff
bison_exe = find_program('bison')
# By convention, all our Bison source files define the api.location.file variable
# so that the generated location file follows the @BASENAME@Location.hh pattern
bison_gen = generator(bison_exe, output : ['@BASENAME@.cc', '@BASENAME@.hh', '@BASENAME@Location.hh'],
                      arguments : ['-W', '-o', '@OUTPUT0@', '@INPUT@'])
bison_src = bison_gen.process('src/macro/Parser.yy', 'src/DynareBison.yy')

incdir = include_directories('src', 'src/macro')

main_src = [ 'src/ComputingTasks.cc',
	     'src/EquationTags.cc',
	     'src/ModelTree.cc',
	     'src/StaticModel.cc',
	     'src/DynamicModel.cc',
	     'src/NumericalConstants.cc',
	     'src/NumericalInitialization.cc',
	     'src/Shocks.cc',
	     'src/SymbolTable.cc',
	     'src/SymbolList.cc',
	     'src/ParsingDriver.cc',
	     'src/DataTree.cc',
	     'src/ModFile.cc',
	     'src/ConfigFile.cc',
	     'src/Statement.cc',
	     'src/ExprNode.cc',
	     'src/VariableDependencyGraph.cc',
	     'src/DynareMain.cc',
	     'src/MacroExpandModFile.cc',
	     'src/Bytecode.cc',
	     'src/ExternalFunctionsTable.cc',
	     'src/ModelEquationBlock.cc',
	     'src/WarningConsolidation.cc',
	     'src/SubModel.cc',
             'src/macro/Driver.cc',
	     'src/macro/Environment.cc',
	     'src/macro/Expressions.cc',
	     'src/macro/Directives.cc' ]

executable('dynare-preprocessor', main_src, flex_src, flexlexer_h, bison_src,
           include_directories : incdir, dependencies : boost_dep,
           link_args : get_option('prefer_static') ? '-static' : '',
           install : true)

## LaTeX stuff

latexmk = find_program('latexmk')

# We have to set TEXINPUTS because the current directory is not the source
# directory when latexmk is invoked (and using the -cd option in combination
# with -outdir/-auxdir does not work in all cases because @OUTDIR@ and
# @PRIVATE_DIR@ can be relative paths)

# Without the -g flag, latexmk remembers a previous build failure and will refuse
# to recompile even if the error has been fixed in the TeX source

# The \graphicspath{} command does not compute directories relative to TEXINPUTS,
# so add these manually
logos_dir = meson.current_source_dir() + '/doc/logos'

custom_target('macroprocessor.pdf',
              output : 'macroprocessor.pdf',
              input : 'doc/macroprocessor/macroprocessor.tex',
              command : [ latexmk, '-pdf', '-g', '-auxdir=@PRIVATE_DIR@', '@INPUT@'],
              env : { 'TEXINPUTS': meson.current_source_dir() + '/doc/macroprocessor:' + logos_dir + ':' },
              install : true,
              install_dir : 'share/doc/dynare')

custom_target('preprocessor.pdf',
              output : 'preprocessor.pdf',
              input : 'doc/preprocessor/preprocessor.tex',
              command : [ latexmk, '-pdf', '-g', '-auxdir=@PRIVATE_DIR@', '@INPUT@'],
              env : { 'TEXINPUTS': meson.current_source_dir() + '/doc/preprocessor:' + logos_dir + ':' },
              install : true,
              install_dir : 'share/doc/dynare')
