ROOT=../..
OBJS=ftsm_var.o

COMMON=..
INCLUDES= -I $(COMMON) -I $(ROOT)/modules


#LOCALMACROS=-include source.defs
#LOCALDEPS=source.defs
#LIBDEP=$(COMMON)/sm_common.o

include $(ROOT)/Makefile.inc

ftsm_var.o: sm_common_vars

sm_common_vars:
	$(MAKE) -C $(COMMON) sm_var.o sm_config.o 

clean: 
	rm -f *.mod $(CHARMM_SOURCES) $(SOURCES) $(OBJS) $(EXE)


.PNONY : sm_common_vars
