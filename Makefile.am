SUBDIRS = macro

BUILT_SOURCES = DynareBison.hh stack.hh position.hh location.hh DynareBison.cc DynareFlex.cc

matlabdir = $(libdir)/matlab

matlab_PROGRAMS = dynare_m

dynare_m_SOURCES = \
	DynareFlex.ll \
	DynareBison.yy \
	ComputingTasks.cc \
	ComputingTasks.hh \
	ModelTree.cc \
	ModelTree.hh \
	StaticModel.cc \
	StaticModel.hh \
	StaticDllModel.cc \
	StaticDllModel.hh \
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
	Statement.cc \
	Statement.hh \
	ExprNode.cc \
	ExprNode.hh \
	MinimumFeedbackSet.cc \
	MinimumFeedbackSet.hh \
	IncidenceMatrix.cc \
	IncidenceMatrix.hh \
	BlockTriangular.cc \
	BlockTriangular.hh \
	DynareMain.cc \
	DynareMain2.cc \
	CodeInterpreter.hh

dynare_m_CPPFLAGS = $(BOOST_CPPFLAGS)
dynare_m_LDFLAGS = $(BOOST_LDFLAGS)
dynare_m_LDADD = macro/libmacro.a

DynareFlex.cc: DynareFlex.ll
	$(LEX) -oDynareFlex.cc DynareFlex.ll

DynareBison.cc DynareBison.hh location.hh stack.hh position.hh: DynareBison.yy
	$(YACC) -o DynareBison.cc DynareBison.yy

all-local:
	cp $(PROGRAMS) ../matlab/

if HAVE_DOXYGEN
html-local:
	$(DOXYGEN)
endif

clean-local:
	cd ../matlab && rm -f $(PROGRAMS)
	rm -rf doc/html/
