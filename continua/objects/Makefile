ROOT=../..
COMMON=..
OBJS=object.o
NAME=continua/objects
DEPS=parser output
SUBDIRS=molecule

include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

object.F90 : function.def subroutine.def 

clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done
