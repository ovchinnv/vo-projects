ROOT=..
OBJS=rng.o

INCLUDES=-I $(ROOT)/modules

include $(ROOT)/Makefile.inc

default: all modules
rng.o : rng.F90 clcg.o constants output
clcg.o : clcg.f clcg.h

output:
	$(MAKE) -C $(ROOT)/output all modules
constants:
	$(MAKE) -C $(ROOT)/constants all modules

.PHONY : constants output

