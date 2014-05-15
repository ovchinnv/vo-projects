# THIS MAKEFILE IS GENERATED AUTOMATICALLY FROM A MASTER SOURCE CODE TREE
ROOT=..
COMMON=$(ROOT)/scommon
OBJS=multicom.o multicom_aux.o ifstack.o
LOCALMACROS=-include source.defs
LOCALDEPS=source.defs
CHARMM_SOURCES=$(OBJS:%.o=%.src)
include $(ROOT)/Makefile.inc
#### RULES:####
charmm : $(CHARMM_SOURCES)


clean ::
	rm -f $(CHARMM_SOURCES) *tmp
