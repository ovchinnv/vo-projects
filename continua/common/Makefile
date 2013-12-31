ROOT=../..
COMMON=..
OBJS=gridsize.o state.o grid.o continuum.o
NAME=continua/common
DEPS=output parser constants files datatypes multicom continua/io continua/objects

include $(COMMON)/Makefile.inc
include $(ROOT)/Makefile.inc

continuum.o: gridsize.o grid.o
state.o: gridsize.o
grid.o: gridsize.o state.o

grid.F90 : gridopts.def

test:: force_look
	$(MAKE) -C test

.SILENT : test
