NAME=tempering
ROOT=../
HEADS=tempering.h
OBJS=tempering_plugin_acemd.o

LOCALDEFS+=-D__BINDC

include $(ROOT)/Makefile.inc

LOCALMACROS+=-include source.defs
LOCALDEPS+=source.defs

include MakeACEMD.inc

tempering_plugin_acemd.o : tempering_plugin_acemd.c aceplug.h $(LOCALDEPS)

#	$(CC) $(CFLAGS) $(LOCALDEFS) $(LOCALMACROS) $< -o $@

clean::
	rm -f $(OBJS) tempering.h
