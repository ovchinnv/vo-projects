ROOT=..
SOURCES=multicom.F90 multicom_aux.F90
OBJS=multicom.o multicom_aux.o
CHARMM_SOURCES=multicom.src multicom_aux.src

LOCALMACROS=-include source.defs
LOCALDEPS=source.defs

INCLUDES=-I $(ROOT)/modules

include $(ROOT)/Makefile.inc

multicom.o:multicom_aux.o output ivector parser parselist
output:
	$(MAKE) -C $(ROOT)/output all modules
ivector:
	$(MAKE) -C $(ROOT)/vectors all modules
parser:
	$(MAKE) -C $(ROOT)/parser all modules
parselist:
	$(MAKE) -C $(ROOT)/parselist all modules

charmm: $(CHARMM_SOURCES)

.PHONY : parser ivector output
