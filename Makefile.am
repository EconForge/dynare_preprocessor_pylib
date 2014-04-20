SUBDIRS = macro

BUILT_SOURCES = DynareBison.hh stack.hh position.hh location.hh DynareBison.cc DynareFlex.cc FlexLexer.h

matlabdir = $(libdir)/matlab

matlab_PROGRAMS = dynare_m

# We don't put BUILT_SOURCES in dynare_m_SOURCES, otherwise DynareBison.o and DynareFlex.o will be linked two times (Automake translates DynareFlex.ll and DynareBison.yy into their respective .o); so BUILT_SOURCES is in EXTRA_DIST
dynare_m_SOURCES = \
	DynareFlex.ll \
	DynareBison.yy \
	ComputingTasks.cc \
	ComputingTasks.hh \
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
	MinimumFeedbackSet.cc \
	MinimumFeedbackSet.hh \
	DynareMain.cc \
	DynareMain2.cc \
	CodeInterpreter.hh \
	ExternalFunctionsTable.cc \
	ExternalFunctionsTable.hh \
	SteadyStateModel.hh \
	SteadyStateModel.cc \
	WarningConsolidation.hh \
	WarningConsolidation.cc \
	ExtendedPreprocessorTypes.hh

# The -I. is for <FlexLexer.h>
dynare_m_CPPFLAGS = $(BOOST_CPPFLAGS) -I.
dynare_m_LDFLAGS = $(BOOST_LDFLAGS)
dynare_m_LDADD = macro/libmacro.a

DynareFlex.cc FlexLexer.h: DynareFlex.ll
	$(LEX) -oDynareFlex.cc DynareFlex.ll
	cp $(LEXINC)/FlexLexer.h .

DynareBison.cc DynareBison.hh location.hh stack.hh position.hh: DynareBison.yy
	$(YACC) -o DynareBison.cc DynareBison.yy

all-local:
	cd ../matlab && $(LN_S) -f $(abs_srcdir)/$(PROGRAMS) $(PROGRAMS)

if HAVE_DOXYGEN
html-local:
	$(DOXYGEN)
endif

clean-local:
	cd ../matlab && rm -f $(PROGRAMS)
	rm -rf doc/html/

EXTRA_DIST = $(BUILT_SOURCES) Doxyfile
