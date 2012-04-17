ROOT=..
OBJS=sm_var.o sm_config.o splines.o sm_util.o

#LOCALMACROS=-include source.defs
#LOCALDEPS=source.defs

include $(ROOT)/Makefile.inc

default: all
all: $(OBJS)
$(OBJS):$(SOURCES)
charmm: $(CHARMM_SOURCES)

clean: 
	rm -f *.mod $(CHARMM_SOURCES) $(SOURCES) $(OBJS) $(EXE)
