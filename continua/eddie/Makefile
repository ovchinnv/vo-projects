ROOT=../..
COMMON=..
OBJS=eddie_grid.o
NAME=continua/eddie
DEPS=output parser constants continua/common

include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

eddie_grid.o :

test:: force_look
	$(MAKE) -C test

.SILENT : test
