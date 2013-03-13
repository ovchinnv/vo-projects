ROOT=../..
COMMON=..
OBJS=gridsize.o state.o grid.o fieldbcdef.o datatypes.o
NAME=continua/common
DEPS=output parser constants files multicom continua/io continua/objects

LOCALDEPS=

include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

state.o: gridsize.o fieldbcdef.o datatypes.o
grid.o: gridsize.o state.o
fieldbcdef.F90 : $(COMMON)/bcutil/bc_util.src
fieldbcdef.o : gridsize.o

clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done

test:: force_look
	$(MAKE) -C test

.SILENT : test
