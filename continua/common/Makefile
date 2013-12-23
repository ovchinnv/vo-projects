ROOT=../..
COMMON=..
OBJS=gridsize.o state.o grid.o fieldbcdef.o datatypes.o continuum.o
NAME=continua/common
DEPS=output parser constants files multicom continua/io continua/objects

include $(COMMON)/Makefile.inc
include $(ROOT)/Makefile.inc

continuum.o: gridsize.o grid.o
state.o: gridsize.o fieldbcdef.o datatypes.o
grid.o: gridsize.o state.o
fieldbcdef.F90 : $(COMMON)/bcutil/bc_util.src
fieldbcdef.o : gridsize.o

grid.F90 : gridopts.def

test:: force_look
	$(MAKE) -C test

.SILENT : test
