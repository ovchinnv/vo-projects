ROOT=..
SOURCES=parser.F90
OBJS=parser.o

INCLUDES=-I $(ROOT)/modules

include $(ROOT)/Makefile.inc

default: all modules
all: $(OBJS)
parser.o : output
output: force_look
	$(MAKE) -C $(ROOT)/output all modules
#$(OBJS):$(SOURCES)
source: $(SOURCES)
clean: 
	rm -f *.mod $(SOURCES) $(OBJS)
