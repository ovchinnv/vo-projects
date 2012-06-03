ROOT=..
OBJS=rng.o
NAME=random

INCLUDES=-I $(ROOT)/modules

include $(ROOT)/Makefile.inc

rng.o : clcg.o constants.stamp output.stamp parser.stamp
clcg.o : clcg.f clcg.h

clean: 
	rm -f *.mod $(SOURCES) $(OBJS) $(ROOT)/$(NAME).stamp
