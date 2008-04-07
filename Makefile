include Makefile.include

ifeq ($(shell uname -o), Cygwin)
	DYNARE_M = dynare_m.exe
	DYNARE_S = dynare_s.exe
else 
	DYNARE_M = dynare_m
	DYNARE_S = dynare_s
endif

ifeq ($(CROSS_WIN32), yes)
	DYNARE_M = dynare_m.exe
	DYNARE_S = dynare_s.exe
endif

COMMON_OBJ = \
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
	Statement.o \
	ExprNode.o \
	ModelNormalization.o \
	ModelBlocks.o \
	BlockTriangular.o \
	Model_Graph.o \
	SymbolGaussElim.o \
	DynareMain.o \
	DynareMain2.o

MATLAB_OBJ = InterfaceMatlab.o

SCILAB_OBJ = InterfaceScilab.o


# Build rules

all: all-recursive $(DYNARE_M) $(DYNARE_S)

all-recursive:
	make -C macro

$(DYNARE_M): $(COMMON_OBJ) $(MATLAB_OBJ) macro/libmacro.a
	$(CPP) $(CPPFLAGS) -o $(DYNARE_M) $(COMMON_OBJ) $(MATLAB_OBJ) -Lmacro -lmacro
	cp $(DYNARE_M) ../matlab/

$(DYNARE_S): $(COMMON_OBJ) $(SCILAB_OBJ) macro/libmacro.a
	$(CPP) $(CPPFLAGS) -o $(DYNARE_S) $(COMMON_OBJ) $(SCILAB_OBJ) -Lmacro -lmacro
	cp $(DYNARE_S) ../scilab/


# Dependencies

-include $(COMMON_OBJ:.o=.P) $(MATLAB_OBJ:.o=.P) $(SCILAB_OBJ:.o=.P)

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
		$(DYNARE_M) \
		$(DYNARE_S)

clean-recursive:
	make -C macro clean

.PHONY: all all-recursive clean clean-recursive
