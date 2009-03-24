CXXFLAGS = -Wall

ifeq ($(shell uname -o), Cygwin)
# Detection of uninitialized variables is buggy in Cygwin and generates spurious warnings
CXXFLAGS += -Wno-uninitialized
CXXFLAGS += -mno-cygwin
endif

ifeq ($(CROSS_WIN32), yes)
CXX = i586-mingw32msvc-g++
AR = i586-mingw32msvc-ar
# Detection of uninitialized variables is buggy in MinGW and generates spurious warnings
CXXFLAGS += -Wno-uninitialized
endif

ifeq ($(DEBUG),yes)
CXXFLAGS += -ggdb
else
CXXFLAGS += -O3
endif

ifeq ($(VALGRIND), yes)
CXXFLAGS = -Wall -O -g -fno-inline
endif

ifeq ($(shell uname -o), Cygwin)
DYNARE_M = dynare_m.exe
else
DYNARE_M = dynare_m
endif

ifeq ($(CROSS_WIN32), yes)
DYNARE_M = dynare_m.exe
endif

MAIN_OBJS = \
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
	IncidenceMatrix.o \
	BlockTriangular.o \
	ModelGraph.o \
	DynareMain.o \
	DynareMain2.o

MACRO_OBJS = \
	macro/MacroFlex.o \
	macro/MacroBison.o \
	macro/MacroDriver.o \
	macro/MacroValue.o


# Build rules

.PHONY: all
all: $(DYNARE_M)

$(DYNARE_M): $(MAIN_OBJS) $(MACRO_OBJS)
	$(CXX) $(CXXFLAGS) -o $(DYNARE_M) $(MAIN_OBJS) $(MACRO_OBJS)
	cp $(DYNARE_M) ../matlab/


# Build rules for Flex and Bison files

DynareFlex.cc: DynareFlex.ll
	flex -oDynareFlex.cc DynareFlex.ll

DynareBison.cc DynareBison.hh location.hh stack.hh position.hh: DynareBison.yy
	bison --verbose -o DynareBison.cc DynareBison.yy

macro/MacroFlex.cc: macro/MacroFlex.ll
	cd macro && flex -oMacroFlex.cc MacroFlex.ll

macro/MacroBison.cc macro/MacroBison.hh macro/location.hh macro/stack.hh macro/position.hh: macro/MacroBison.yy
	cd macro && bison --verbose -o MacroBison.cc MacroBison.yy


# Dependencies

# General rule for creating per-source dependencies Makefile
# We use -MG to avoid failing on generated headers (MacroBison.hh, DynareBison.hh)
# As a consequence, these headers are included without path-prefix
%.d: %.cc
	@set -e; rm -f $@; \
	 $(CXX) -MM -MG $(CPPFLAGS) $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	 rm -f $@.$$$$

# As DynareBison.hh, this file is included in the .d files without its path (since it is generated), so we force the path
vpath MacroBison.hh macro

-include $(MAIN_OBJS:.o=.d)
-include $(MACRO_OBJS:.o=.d)


# Clean

.PHONY: clean
clean:
	rm -f *.o *.d *~ \
		DynareFlex.cc \
		DynareBison.output \
		DynareBison.cc \
		position.hh \
		stack.hh \
		location.hh \
		DynareBison.hh \
		$(DYNARE_M)
	cd macro && rm -f *.o *.d *~ \
		MacroFlex.cc \
		MacroBison.output \
		MacroBison.cc \
		MacroBison.hh \
		location.hh \
		stack.hh \
		position.hh
