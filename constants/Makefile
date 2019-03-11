ROOT=..
OBJS=constants.o
MODS=constants.mod
NAME=constants

# make sure not to use ofast in single precision
# this will cause error tolerances to appear zero
LOCALDEFS+=-O3 -fprotect-parens
TESTOBJS=$(OBJS)
include $(ROOT)/Makefile.inc

