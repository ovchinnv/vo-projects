NAME=tempering
ROOT=..
HEADS=tempering.h
OBJS= tempering_plugin.o tempering.o griddata.o energy_grid.o
COBJS=tempering_plugin_acemd.o

DEPS=files output parser random constants vectors

PYTHONH=-I /usr/include/python3.8/
PYTHONINSTALL=/usr/lib/python3.8/site-packages
LPYTHON=-lpython3.8

#PYTHONH=-I /usr/include/python3.6m/
#PYTHONINSTALL=/usr/lib/python3.6/site-packages
#LPYTHON=-lpython3.6m

#PYTHONH=-I /usr/include/python3.2mu/
#PYTHONINSTALL=/usr/lib64/python3.2/site-packages
#LPYTHON=-lpython3.2mu

#PYTHONH=-I /usr/include/python2.7/
#PYTHONINSTALL=/usr/lib64/python2.7/site-packages
#LPYTHON=-lpython2.7

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
aceplug.h:
	wget raw.githubusercontent.com/Acellera/plumed1/master/plumed-1.3/ACEMD/aceplug.h
tempering_plugin.o : tempering.o
tempering.o : energy_grid.o
energy_grid.o : griddata.o

tempering.F90 : parser/parser_macro2.src
