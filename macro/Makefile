include ../Makefile.include

OBJ = MacroFlex.o MacroBison.o MacroDriver.o

libmacro.a: $(OBJ)
	$(AR) crs libmacro.a $(OBJ)

-include $(OBJ:.o=.P)

MacroFlex.cc: MacroFlex.ll MacroBison.hh MacroDriver.hh
	flex -oMacroFlex.cc MacroFlex.ll

MacroBison.cc MacroBison.hh: MacroBison.yy MacroDriver.hh
	bison --verbose -o MacroBison.cc MacroBison.yy

clean: 
	rm -f *.o *.P *~ \
		MacroFlex.cc \
		MacroBison.output \
		MacroBison.cc \
		MacroBison.hh \
		location.hh \
		stack.hh \
		position.hh \
		libmacro.a

.PHONY: clean
