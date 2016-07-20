NAME=tempering
ROOT=../
HEADS=tempering.h
OBJS=tempering_plugin_acemd.o tempering_acemd.o tempering.o griddata.o energy_grid.o
DEPS=files output parser random

LOCALDEFS+=-D__BINDC

include $(ROOT)/Makefile.inc

LOCALMACROS+=-include source.defs
LOCALDEPS+=source.defs

include MakeACEMD.inc

$(SOURCES) : $(LOCALDEPS)

tempering_plugin_acemd.o : tempering_plugin_acemd.c aceplug.h
tempering_acemd.o : tempering.o
tempering.o : energy_grid.o
energy_grid.o : griddata.o

tempering.F90 : parser/parser_macro2.src
