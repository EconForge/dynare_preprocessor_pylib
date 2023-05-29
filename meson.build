# TODO for reaching feature parity with current autotools-based build system:
# - Detection of Boost Graph library and corresponding include directory
# - Detection of FlexLexer.h and corresponding include directory
# - Build and install documentation (with a configuration option to disable it explicitly)
# - Link with pthread with old glibc
# - Add -Wold-style-cast flag when building flex-generated files

project('dynare-preprocessor', 'cpp',
        version : '6-unstable',
        default_options : [ 'cpp_std=gnu++20', 'warning_level=2' ])

add_global_arguments('-DPACKAGE_VERSION="' + meson.project_version() + '"', language : 'cpp')

flex_exe = find_program('flex')
# The -Ca flag comes from hitting a hard-coded size limit.
# Partial explanation: https://www.owlfolio.org/possibly-useful/flex-input-scanner-rules-are-too-complicated
# There is a Debian bug report about this: https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=642040
flex_gen = generator(flex_exe, output : '@BASENAME@.cc',
                     arguments : ['-Ca', '-o', '@OUTPUT@', '@INPUT@'])

bison_exe = find_program('bison')
# By convention, all our Bison source files define the api.location.file variable
# so that the generated location file follows the @BASENAME@Location.hh pattern
bison_gen = generator(bison_exe, output : ['@BASENAME@.cc', '@BASENAME@.hh', '@BASENAME@Location.hh'],
                      arguments : ['-W', '-o', '@OUTPUT0@', '@INPUT@'])

flex_src = flex_gen.process('src/macro/Tokenizer.ll', 'src/DynareFlex.ll')
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

executable('dynare-preprocessor', main_src, flex_src, bison_src,
           include_directories : incdir, install : true)
