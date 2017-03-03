NAME=watershell
ROOT=..
HEADS=watershell.h
OBJS= watershell.o watershell_plugin.o watershell_plugin_acemd.o

DEPS=files dynamol continua output parser constants erf

LIB=libwatershell
LIBS+=$(LIB).so

include linkobjs.inc
include $(ROOT)/Makefile.inc

LOCALMACROS+=-include source.defs -D__BINDC -D__erf=erfo7
LOCALDEPS+=source.defs

default : plugin
plugin : $(LIB).so
$(LIB).so : $(OBJS) $(COBJS) $(LOBJS) Makefile.inc

$(SOURCES) : $(LOCALDEPS)

watershell_plugin.o : watershell.o
watershell.o :

watershell.F90 watershell_plugin.F90 : parser/parser_macro2.src
