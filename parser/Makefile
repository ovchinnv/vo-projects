ROOT=..
OBJS=parser.o

INCLUDES=-I $(ROOT)/modules

include $(ROOT)/Makefile.inc

default: all modules
all: $(OBJS)
parser.o : output
output:
	$(MAKE) -C $(ROOT)/output all modules

.PHONY : output