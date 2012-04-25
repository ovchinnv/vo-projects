ROOT=..
OBJS=parselist.o

INCLUDES=-I $(ROOT)/modules

include $(ROOT)/Makefile.inc

default: all modules
all: $(OBJS)
iparselist.o : output ivector
ivector:
	$(MAKE) -C $(ROOT)/vectors all modules
$(OBJS):$(SOURCES)
source: $(SOURCES)
charmm: $(CHARMM_SOURCES)
clean: 
	rm -f *.mod $(SOURCES) $(CHARMM_SOURCES)  $(OBJS)

.PHONY : output
