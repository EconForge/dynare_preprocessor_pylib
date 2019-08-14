noinst_LIBRARIES = libmacro.a

BUILT_SOURCES = Parser.hh location.hh Parser.cc Tokenizer.cc

# We don't put BUILT_SOURCES in libmacro_a_SOURCES, otherwise Parser.o and Tokenizer.o will be linked two times (Automake translates Tokenizer.ll and Parser.yy into their respective .o); so BUILT_SOURCES is in EXTRA_DIST
libmacro_a_SOURCES = \
	Tokenizer.ll \
	Parser.yy \
	ForwardDeclarationsAndEnums.hh \
	Driver.cc \
	Driver.hh \
	Environment.cc \
	Environment.hh \
	Expressions.cc \
	Expressions.hh \
	Directives.cc \
	Directives.hh

EXTRA_DIST = $(BUILT_SOURCES)

# The -I.. is for <FlexLexer.h>
libmacro_a_CPPFLAGS = $(BOOST_CPPFLAGS) -I..

Tokenizer.cc: Tokenizer.ll
	$(LEX) -o Tokenizer.cc Tokenizer.ll

libmacro_a-Tokenizer.$(OBJEXT): CXXFLAGS += -Wno-old-style-cast

# We do not use a multiple target rule for Bison, because otherwise it will be
# run several times in parallel builds
Parser.cc: Parser.yy
	$(YACC) -W -o $@ $<

Parser.hh location.hh: Parser.cc
