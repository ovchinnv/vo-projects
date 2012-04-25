ROOT=../..
OBJS=ftsm_var.o ftsm_rex.o 

COMMON=..
INCLUDES= -I $(COMMON) -I $(ROOT)/modules

LOCALMACROS=-include source.defs -include ../source.defs
LOCALDEPS=source.defs ../source.defs

include $(ROOT)/Makefile.inc

ftsm.o: ftsm_var.o ftsm_rex.o
ftsm_var.o: sm_common_vars
ftsm_rex.o: constants ivector ftsm_var.o

sm_common_vars:
	$(MAKE) -C $(COMMON) sm_var.o sm_config.o 

constants:
	$(MAKE) -C $(ROOT)/constants all modules 

ivector:
	$(MAKE) -C $(ROOT)/vectors all modules

charmm: $(CHARMM_SOURCES)


.PHONY : sm_common_vars constants ivector
