noinst_LIBRARIES = libmacro.a

BUILT_SOURCES = MacroBison.hh stack.hh position.hh location.hh MacroBison.cc MacroFlex.cc

# We don't put BUILT_SOURCES in libmacro_a_SOURCES, otherwise MacroBison.o and MacroFlex.o will be linked two times (Automake translates MacroFlex.ll and MacroBison.yy into their respective .o); so BUILT_SOURCES is in EXTRA_DIST
libmacro_a_SOURCES = \
	MacroFlex.ll \
	MacroBison.yy \
	MacroDriver.cc \
	MacroDriver.hh \
	MacroValue.cc \
	MacroValue.hh

EXTRA_DIST = $(BUILT_SOURCES)

# The -I.. is for <FlexLexer.h>
libmacro_a_CPPFLAGS = -I..

MacroFlex.cc: MacroFlex.ll
	$(LEX) -oMacroFlex.cc MacroFlex.ll

MacroBison.cc MacroBison.hh location.hh stack.hh position.hh: MacroBison.yy
	$(YACC) -o MacroBison.cc MacroBison.yy
