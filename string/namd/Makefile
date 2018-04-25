NAME=string/namd
ROOT=../..
COMMON=..
HEADS=GlobalMasterSM.h

LOCALDEFS+=-D__BINDC -D__CPLUSPLUS

include $(ROOT)/Makefile.inc

LOCALMACROS+=-include source.defs
LOCALMACROS+=-include $(COMMON)/source.defs
LOCALDEPS+=source.defs
LOCALDEPS+=$(COMMON)/source.defs

include $(COMMON)/Makefile.inc

GlobalMasterSM.h:GlobalMasterSM.H ../plugin/struna.H $(LOCALDEPS) $(GLOBALDEPS)
	$(FPP) $(FPPFLAGS) $(DEFS) $(LOCALDEFS) $(GLOBALMACROS) $(LOCALMACROS) -P $< > $@
	cat source.defs >> GlobalMasterSM.h

clean::
	rm -f $(OBJS) $(HEADS)
