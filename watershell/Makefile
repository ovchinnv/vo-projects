NAME=watershell
ROOT=..
HEADS=watershell.h
OBJS= watershell.o
COBJS=watershell_plugin_acemd.o

DEPS=files dynamol continua output parser constants

include $(ROOT)/Makefile.inc

LOCALMACROS+=-include source.defs
LOCALDEPS+=source.defs

$(SOURCES) : $(LOCALDEPS)

watershell.o :

watershell.F90 : parser/parser_macro2.src
