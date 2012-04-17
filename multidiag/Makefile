ROOT=..
OBJS=multidiag.o

include $(ROOT)/Makefile.inc

default: all modules
all: $(OBJS)
$(OBJS):$(SOURCES)
source: $(SOURCES)
charmm: $(CHARMM_SOURCES)
clean: 
	rm -f *.mod $(SOURCES) $(CHARMM_SOURCES) $(OBJS)
