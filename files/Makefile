ROOT=..
OBJS=files.o

INCLUDES=-I $(ROOT)/modules

include $(ROOT)/Makefile.inc

default: all modules
files.o : parser ivector output
output:
	$(MAKE) -C $(ROOT)/output all modules
ivector:
	$(MAKE) -C $(ROOT)/vectors all modules
parser:
	$(MAKE) -C $(ROOT)/parser all modules

.PHONY : parser ivector output

