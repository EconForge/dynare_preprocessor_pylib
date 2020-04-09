SUBDIRS = macro

BUILT_SOURCES = DynareBison.hh location.hh DynareBison.cc DynareFlex.cc FlexLexer.h

bin_PROGRAMS = dynare_m

# We don't put BUILT_SOURCES in dynare_m_SOURCES, otherwise DynareBison.o and DynareFlex.o will be linked two times (Automake translates DynareFlex.ll and DynareBison.yy into their respective .o); so BUILT_SOURCES is in EXTRA_DIST
dynare_m_SOURCES = \
	DynareFlex.ll \
	DynareBison.yy \
	ComputingTasks.cc \
	ComputingTasks.hh \
	EquationTags.cc \
	EquationTags.hh \
	ModelTree.cc \
	ModelTree.hh \
	StaticModel.cc \
	StaticModel.hh \
	DynamicModel.cc \
	DynamicModel.hh \
	NumericalConstants.cc \
	NumericalConstants.hh \
	NumericalInitialization.cc \
	NumericalInitialization.hh \
	Shocks.cc \
	Shocks.hh \
	SigmaeInitialization.cc \
	SigmaeInitialization.hh \
	SymbolTable.cc \
	SymbolTable.hh \
	SymbolList.cc \
	SymbolList.hh \
	ParsingDriver.cc \
	ParsingDriver.hh \
	DataTree.cc \
	DataTree.hh \
	ModFile.cc \
	ModFile.hh \
	ConfigFile.cc \
	ConfigFile.hh \
	Statement.cc \
	Statement.hh \
	ExprNode.cc \
	ExprNode.hh \
	VariableDependencyGraph.cc \
	VariableDependencyGraph.hh \
	DynareMain.cc \
	MacroExpandModFile.cc \
	CodeInterpreter.hh \
	ExternalFunctionsTable.cc \
	ExternalFunctionsTable.hh \
	ModelEquationBlock.hh \
	ModelEquationBlock.cc \
	WarningConsolidation.hh \
	WarningConsolidation.cc \
	ExtendedPreprocessorTypes.hh \
	SubModel.cc \
	SubModel.hh


ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = \
	Doxyfile \
	$(BUILT_SOURCES)

# The -I. is for <FlexLexer.h>
dynare_m_CPPFLAGS = $(BOOST_CPPFLAGS) -I.
dynare_m_LDFLAGS = $(AM_LDFLAGS) $(BOOST_LDFLAGS)
dynare_m_LDADD = macro/libmacro.a -lstdc++fs

# -Ca flag comes from hitting a hard-coded size limit.
# Partial explanation: https://www.owlfolio.org/possibly-useful/flex-input-scanner-rules-are-too-complicated
# There is a Debian bug report about this: https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=642040
DynareFlex.cc: DynareFlex.ll
	$(LEX) -Ca -o DynareFlex.cc DynareFlex.ll

FlexLexer.h:
	cp $(LEXINC)/FlexLexer.h . || test -f ./FlexLexer.h

dynare_m-DynareFlex.$(OBJEXT): CXXFLAGS += -Wno-old-style-cast

# We do not use a multiple target rule for Bison, because otherwise it will be
# run several times in parallel builds
DynareBison.cc: DynareBison.yy
	$(YACC) -W -o $@ $<

DynareBison.hh location.hh: DynareBison.cc

all-local: $(PROGRAMS)

if HAVE_DOXYGEN
html-local:
	$(DOXYGEN)
endif

clean-local:
	rm -rf doc/html/
