ROOT=../../..
COMMON=../..
OBJS=PDB.o molecule.o
DEPS=dynamol output parser erf
NAME=continua/objects/molecule
#
include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc
#
molecule.o: PDB.o
molecule.F90: mol_macro1.def
