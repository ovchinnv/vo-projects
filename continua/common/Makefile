ROOT=../..
COMMON=..
OBJS=gridsize.o continuum.o grid.o
NAME=continua/common
DEPS=output parser constants files state multicom continua/io continua/objects

LOCALDEFS+=-D__DEBUG

include $(COMMON)/Makefile.inc
include $(ROOT)/Makefile.inc

continuum.o: gridsize.o grid.o
grid.o: gridsize.o
grid.F90 : gridopts.def $(ROOT)/parser/parser_macro1.src

test:: force_look
	$(MAKE) -C test

.SILENT : test
