SHELL=/bin/sh
GLOBALMACROS=-include $(ROOT)/source.defs
MSG=-include $(ROOT)/source.msg
CHMSG=-include $(ROOT)/charmm.msg
FPP=cpp
FPPFLAGS = -P -C -I $(ROOT) -nostdinc # -nostdinc : avoid including files automatically
CHARMM_MAJOR_VERSION=38
#### RULES:####
%.src:%.ftn $(LOCALDEPS) $(GLOBALDEPS) $(ROOT)/charmm.msg
	$(FPP) $(FPPFLAGS)         -D__CHARMM -D__CHARMM_VER=$(CHARMM_MAJOR_VERSION) $(GLOBALMACROS) $(LOCALMACROS) $(CHMSG) -P $< > $@

