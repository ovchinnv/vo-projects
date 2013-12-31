ROOT=../
OBJS=fieldbcdef.o datatypes.o
NAME=datatypes
DEPS=parser

include $(ROOT)/Makefile.inc

datatypes.o: fieldbcdef.o
fieldbcdef.F90 : $(ROOT)/datatypes/bcutil/bc_util.src

test:: force_look
	$(MAKE) -C test

.SILENT : test
