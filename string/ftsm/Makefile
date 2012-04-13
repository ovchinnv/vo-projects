ROOT=../..
CHARMM_SOURCES=ftsm_var.src
include $(ROOT)/Makefile.inc

COMMON=..

#LOCALMACROS=-include source.defs
LOCALDEPS=source.defs
$(COMMON)/sm_common.o

OBJS=ftsm_var.o
SOURCES=ftsm_var.F90

default: all $(GLOBALDEPS) $(LOCALDEPS)
all: $(OBJS)
%.F90:%.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(MSG) -P $< > $@
%.src:%.ftn
	$(FPP) $(FPPFLAGS)         -D__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(MSG) -P $< > $@

ftsm_var.o   : ftsm_var.F90
ftsm_var.F90 : ftsm_var.ftn

charmm: $(CHARMM_SOURCES) $(GLOBALDEPS) $(LOCALDEPS)
ftsm_var.src: ftsm_var.ftn

clean: 
	rm -f *.mod $(CHARMM_SOURCES) $(SOURCES) $(OBJS) $(EXE)
