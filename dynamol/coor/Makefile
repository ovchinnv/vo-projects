ROOT=../..
COMMON=..
OBJS = mol_formats.o charmmio.o freeio.o pdbio.o corio.o
NAME=dynamol/coor
DEPS=parser output dynamol/struc

LOCALMACROS=-include $(COMMON)/source.defs -include $(ROOT)/source.defs 
LOCALDEPS=$(COMMON)/source.defs $(ROOT)/source.defs 

include $(ROOT)/Makefile.inc
corio.o : mol_formats.o charmmio.o freeio.o pdbio.o
