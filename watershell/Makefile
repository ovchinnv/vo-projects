NAME=watershell
ROOT=..
HEADS=watershell.h
OBJS= watershell.o watershell_plugin.o watershell_plugin_acemd.o ugrid.o

DEPS=files dynamol continua output parser constants erf lu unitcell

LIB=libwatershell.so # libwatershell.a
LIBS+=$(LIB)

include linkobjs.inc
include $(ROOT)/Makefile.inc

LOCALMACROS+=-include source.defs -D__BINDC -include erf_as.def #-D__erf=erfo7
LOCALDEPS+=source.defs erf_as.def

default : plugin
plugin : $(LIB)
$(LIB) : $(OBJS) $(COBJS) $(LOBJS) linkobjs.inc

$(SOURCES) : $(LOCALDEPS)

watershell.o : ugrid.o

watershell_plugin.o : watershell.o

watershell.F90 watershell_plugin.F90 : parser/parser_macro2.src parser/parser_macro3.src watershell_rmft_aux1.ftn

watershell_plugin_acemd.o : aceplug.h

aceplug.h:
	wget raw.githubusercontent.com/Acellera/plumed1/master/plumed-1.3/ACEMD/aceplug.h
