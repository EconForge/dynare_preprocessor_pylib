noinst_LIBRARIES = libmacro.a

BUILT_SOURCES = MacroBison.hh stack.hh position.hh location.hh MacroBison.cc MacroFlex.cc

libmacro_a_SOURCES = \
	MacroFlex.ll \
	MacroBison.yy \
	MacroDriver.cc \
	MacroDriver.hh \
	MacroValue.cc \
	MacroValue.hh

MacroFlex.cc: MacroFlex.ll
	$(LEX) -oMacroFlex.cc MacroFlex.ll

MacroBison.cc MacroBison.hh location.hh stack.hh position.hh: MacroBison.yy
	$(YACC) -o MacroBison.cc MacroBison.yy
