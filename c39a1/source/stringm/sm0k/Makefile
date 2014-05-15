# THIS MAKEFILE IS GENERATED AUTOMATICALLY FROM A MASTER SOURCE CODE TREE
ROOT=..
COMMON=$(ROOT)/scommon
OBJS=sm0k.o
LOCALMACROS=-include $(COMMON)/source.defs
LOCALDEPS=$(ROOT)/mpitype.def $(COMMON)/source.defs 
CHARMM_SOURCES=$(OBJS:%.o=%.src)
include $(ROOT)/Makefile.inc
#### RULES:####
charmm : $(CHARMM_SOURCES)


clean ::
	rm -f $(CHARMM_SOURCES) *tmp
