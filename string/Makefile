ROOT=..
OBJS=sm_var.o sm_config.o splines.o sm_util.o sm_main.o

LOCALMACROS=-include source.defs
LOCALDEPS=mpitype.def source.defs

include $(ROOT)/Makefile.inc

default: all
all: $(OBJS)
$(OBJS):$(SOURCES)
charmm: $(CHARMM_SOURCES)
source: $(SOURCES)
clean: 
	rm -f *.mod $(CHARMM_SOURCES) $(SOURCES) $(OBJS) $(EXE)
