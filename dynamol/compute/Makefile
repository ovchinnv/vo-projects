ROOT=../..
OBJS=angles.o bonds.o dihedrals.o
DIRS=../param ../struc
INCLUDES=-I $(ROOT)/modules
include $(ROOT)/Makefile.inc

modules : interface.h
angles.o dihedrals.o : constants $(DIRS)
#angles.o    : anglpar.o constants
#bonds.o     : bondpar.o constants
#dihedrals.o : dihepar.o constants
