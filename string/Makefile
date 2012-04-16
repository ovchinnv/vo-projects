ROOT=..
SOURCES=sm_var.F90 sm_config.F90 splines.F90 sm_util.F90
OBJS=sm_var.o sm_config.o splines.o sm_util.o
CHARMM_SOURCES=sm_var.src sm_config.src splines.src sm_util.src

#LOCALMACROS=-include source.defs
#LOCALDEPS=source.defs

include $(ROOT)/Makefile.inc

default: all
all: $(OBJS)
$(OBJS):$(SOURCES)
charmm: $(CHARMM_SOURCES)

clean: 
	rm -f *.mod $(CHARMM_SOURCES) $(SOURCES) $(OBJS) $(EXE)
