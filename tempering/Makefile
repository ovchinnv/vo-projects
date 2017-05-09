NAME=tempering
ROOT=..
HEADS=tempering.h
OBJS= tempering_plugin.o tempering.o griddata.o energy_grid.o
COBJS=tempering_plugin_acemd.o

DEPS=files output parser random constants vectors

PYTHONH=-I /usr/include/python3.6m/
PYTHONINSTALL=/usr/lib/python3.6/site-packages
LPYTHON=-lpython3.6m
LTCL=-ltcl

SWIGC=swig
SFLAGS=

LOCALDEFS+=-D__BINDC $(PYTHONH)

include $(ROOT)/Makefile.inc

LOCALMACROS+=-include source.defs
LOCALDEPS+=source.defs

include Makefile.inc

$(SOURCES) : $(LOCALDEPS)

tempering_plugin_acemd.o : tempering_plugin_acemd.c aceplug.h
tempering_plugin.o : tempering.o
tempering.o : energy_grid.o
energy_grid.o : griddata.o

tempering.F90 : parser/parser_macro2.src
