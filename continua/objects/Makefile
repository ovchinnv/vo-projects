ROOT=../..
COMMON=..
OBJS=object.o
NAME=continua/objects
DEPS=parser output
SUBDIRS=molecule

include $(COMMON)/Makefile.inc
include $(ROOT)/Makefile.inc

object.F90 : routine.def 

clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done
