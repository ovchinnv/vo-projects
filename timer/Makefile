ROOT=..
OBJS=timer.o

INCLUDES=-I $(ROOT)/modules

include $(ROOT)/Makefile.inc

default: all modules

timer.o: output parser ivector

output:
	$(MAKE) -C $(ROOT)/output all modules
ivector:
	$(MAKE) -C $(ROOT)/vectors all modules
parser:
	$(MAKE) -C $(ROOT)/parser all modules

