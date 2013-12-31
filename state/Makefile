ROOT=../
OBJS=fieldbcdef.o datatypes.o state.o
NAME=state
DEPS=output constants parser

include $(ROOT)/Makefile.inc

datatypes.o: fieldbcdef.o
fieldbcdef.F90 : $(ROOT)/state/bcutil/bc_util.src

test:: force_look
	$(MAKE) -C test

.SILENT : test
