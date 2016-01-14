ROOT=../../..
STRING=../../
SMCV=../
OBJS=plugin-smcv.o minimal_system.o

NAME=string/smcv/acemd

LOCALDEFS+=-D__BINDC

include $(ROOT)/Makefile.inc

LOCALMACROS+=-include $(SMCV)/source.defs -include $(ROOT)/source.defs -include source.defs
LOCALDEPS+=$(SMCV)/source.defs $(ROOT)/source.defs source.defs

plugin-smcv.o : plugin-smcv.c aceplug.h smcv.h $(LOCALDEPS)
	$(CC) $(CFLAGS) $(LOCALDEFS) $(LOCALMACROS) $< -o $@

clean::
	rm -f $(OBJS) smcv.h
