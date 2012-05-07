ROOT=..
OBJS=rng.o

INCLUDES=-I $(ROOT)/modules

include $(ROOT)/Makefile.inc

rng.o : clcg.o constants output parser
clcg.o : clcg.f clcg.h

mod:
	/bin/sh -c 'for file in "rng.mod"; do  ln -fs `pwd`/$$file ${ROOT}/include/$$file; done'
clean: 
	rm -f *.mod $(SOURCES) $(OBJS)
