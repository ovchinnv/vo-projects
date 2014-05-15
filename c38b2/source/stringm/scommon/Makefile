# THIS MAKEFILE IS GENERATED AUTOMATICALLY FROM A MASTER SOURCE CODE TREE
ROOT=..
COMMON=$(ROOT)/scommon
OBJS=sm_var.o sm_config.o splines.o sm_util.o sm_main.o tsp.o
LOCALMACROS=-include source.defs
LOCALDEPS=$(ROOT)/mpitype.def source.defs 
CHARMM_SOURCES=$(OBJS:%.o=%.src)
include $(ROOT)/Makefile.inc
#### RULES:####
charmm : $(CHARMM_SOURCES)


clean ::
	rm -f $(CHARMM_SOURCES) *tmp
