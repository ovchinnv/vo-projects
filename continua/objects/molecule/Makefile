ROOT=../../..
COMMON=../..
OBJS=molecule.o
DEPS=dynamol output parser erf
NAME=continua/objects/molecule

LOCALMACROS+=-include $(ROOT)/erf/erf_as.def -include source.def
LOCALDEPS+=$(ROOT)/erf/erf_as.def source.def

#
include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc
LOCALDEFS+=-fprotect-parens
#
molecule.F90: mol_macro1.src molecule_grid_objects.src source.def
