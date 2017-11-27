ROOT=..
COMMON=.
OBJS=
NAME=continua
SUBDIRS=io common objects chest # eddie
DEPS=state

LOCALMACROS=-include source.defs
LOCALDEPS=$(ROOT)/mpitype.def source.defs 
LOCALDEPS+=outmacro1.def

include $(ROOT)/Makefile.inc
#include $(COMMON)/Makefile.inc

clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done
