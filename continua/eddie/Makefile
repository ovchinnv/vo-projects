ROOT=../..
COMMON=..
OBJS=eddie_grid.o
NAME=continua/eddie
DEPS=output parser constants continua/common

include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

eddie_grid.o :
eddie_grid.F90: $(ROOT)/parser/parser_macro1.src

test:: force_look
	$(MAKE) -C test

.SILENT : test
