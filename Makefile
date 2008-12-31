include Makefile.include

ifeq ($(shell uname -o), Cygwin)
	DYNARE_M = dynare_m.exe
else
	DYNARE_M = dynare_m
endif

ifeq ($(CROSS_WIN32), yes)
	DYNARE_M = dynare_m.exe
endif

OBJS = \
	DynareFlex.o \
	DynareBison.o \
	ComputingTasks.o \
	ModelTree.o \
	NumericalConstants.o \
	NumericalInitialization.o \
	Shocks.o \
	SigmaeInitialization.o \
	SymbolTable.o \
	SymbolList.o \
	VariableTable.o \
	ParsingDriver.o \
	DataTree.o \
	ModFile.o \
	MatlabFile.o \
	Statement.o \
	ExprNode.o \
	ModelNormalization.o \
	ModelBlocks.o \
	IncidenceMatrix.o \
	BlockTriangular.o \
	Model_Graph.o \
	SymbolGaussElim.o \
	DynareMain.o \
	DynareMain2.o


# Build rules

all: all-recursive $(DYNARE_M)

all-recursive:
	make -C macro

$(DYNARE_M): $(OBJS) macro/libmacro.a
	$(CXX) $(CXXFLAGS) -o $(DYNARE_M) $(OBJS) -Lmacro -lmacro
	cp $(DYNARE_M) ../matlab/


# Dependencies

-include $(OBJS:.o=.P)

DynareFlex.cc: DynareFlex.ll include/DynareBison.hh include/ParsingDriver.hh
	flex -oDynareFlex.cc DynareFlex.ll

DynareBison.cc include/DynareBison.hh: DynareBison.yy include/ParsingDriver.hh
	bison --verbose -o DynareBison.cc DynareBison.yy
	mv DynareBison.hh location.hh stack.hh position.hh include/


# Clean

clean: clean-recursive
	rm -f *.o *.P \
		*~ include/*~ \
		DynareFlex.cc \
		DynareBison.output \
		DynareBison.cc \
		include/position.hh \
		include/stack.hh \
		include/location.hh \
		include/DynareBison.hh \
		$(DYNARE_M)

clean-recursive:
	make -C macro clean

.PHONY: all all-recursive clean clean-recursive
