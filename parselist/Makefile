ROOT=..
OBJS=parselist.o

INCLUDES=-I $(ROOT)/modules

include $(ROOT)/Makefile.inc

default: all modules
iparselist.o : ivector
ivector:
	$(MAKE) -C $(ROOT)/vectors all modules
charmm: $(CHARMM_SOURCES)

.PHONY : output
