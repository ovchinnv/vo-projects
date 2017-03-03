NAME=dynamo
ROOT=..
HEADS=plugin.h
OBJS= plugin_master.o plugin_acemd.o

DEPS=files dynamol continua output parser constants

LIB=lib$(NAME).so
LIBS+=$(LIB)

LOCALMACROS+=-include source.defs -D__BINDC 
LOCALDEPS+=source.defs

include plugins.inc
include watershell.inc
include $(ROOT)/Makefile.inc

default : all
all : $(LIB)
$(LIB) : $(OBJS) $(COBJS) $(LOBJS)

$(SOURCES) : $(LOCALDEPS)

plugin_acemd.o : aceplug.h
plugin_master.F90 : parser/parser_macro2.src
