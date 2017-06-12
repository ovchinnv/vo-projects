NAME=plugin_master/namd
ROOT=../..
COMMON=..
HEADS=GlobalMasterDY.h

LOCALDEFS+=-D__BINDC -D__CPLUSPLUS

include $(ROOT)/Makefile.inc

LOCALMACROS+=-include $(COMMON)/source.defs
LOCALDEPS+=$(COMMON)/source.defs

GlobalMasterDY.h:GlobalMasterDY.H ../plugin.H ../source.defs $(LOCALDEPS) $(GLOBALDEPS)
	$(FPP) $(FPPFLAGS) $(DEFS) $(LOCALDEFS) $(GLOBALMACROS) $(LOCALMACROS) -P $< > $@
	cat ../source.defs >> GlobalMasterDY.h

clean::
	rm -f $(OBJS) $(HEADS)
