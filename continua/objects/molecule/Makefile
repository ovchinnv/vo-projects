ROOT=../../..
COMMON=../..
OBJS=molecule.o
DEPS=dynamol output parser erf
NAME=continua/objects/molecule

LOCALMACROS+=-include $(ROOT)/erf/erf_as.def
LOCALDEPS+=$(ROOT)/erf/erf_as.def

#
include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc
#
molecule.F90: mol_macro1.src
