ROOT=../../..
COMMON=../..
OBJS=PDB.o molecule.o
DEPS=dynamol output parser
NAME=chest/objects/molecule
#
include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc
#
molecule.o: PDB.o
molecule.F90: mol_macro1.def
