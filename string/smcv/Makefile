ROOT=../..
OBJS=cv_types.o cv_common.o

COMMON=..
INCLUDES= -I $(COMMON) -I $(ROOT)/modules

LOCALMACROS=-include source.defs -include ../source.defs
LOCALDEPS=source.defs ../source.defs

include $(ROOT)/Makefile.inc

cv_common.o: constants ivector sm_common_vars

sm_common_vars:
	$(MAKE) -C $(COMMON) sm_var.o sm_config.o 

constants:
	$(MAKE) -C $(ROOT)/constants all modules 

ivector:
	$(MAKE) -C $(ROOT)/vectors all modules

multidiag:
	$(MAKE) -C $(ROOT)/multidiag modules

lu:
	$(MAKE) -C $(ROOT)/lu modules

charmm: $(CHARMM_SOURCES)


.PHONY : sm_common_vars constants ivector multidiag lu
