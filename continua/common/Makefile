ROOT=../..
COMMON=..
OBJS=gridsize.o continuum.o
NAME=continua/common
DEPS=output parser constants files state multicom continua/io continua/objects

include $(COMMON)/Makefile.inc
include $(ROOT)/Makefile.inc

continuum.o: gridsize.o grid.o
grid.o: gridsize.o

grid.F90 : gridopts.def

test:: force_look
	$(MAKE) -C test

.SILENT : test
