NAME=griddata
ROOT=..
OBJS=gridscalartype.o gridscalar_log.o gridscalar.o gridscalar_beta.o

DEPS=output constants

LOCALMACROS+=-include source.defs
LOCALDEPS+=source.defs

include $(ROOT)/Makefile.inc
$(SOURCES) : $(LOCALDEPS)

gridscalar_log.o gridscalar.o : gridscalartype.o
gridscalar_log.F90: gridscalartemplate.ftn $(LOCALDEPS) $(GLOBALDEPS) $(HEADS)
	sed 's/\([^\#]\)\(\/\/\)/\1\^\^/g' $< | $(FPP) $(FPPFLAGS) -D'_TAG=_log' \
	-D'__XMAP(_X)=log(_X)' -D'__XIMAP(_X)=exp(_X)' -D'__XGTZERO' $(DEFS) $(LOCALDEFS) -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(MSG) $(INCLUDES) -P | sed 's/\^\^/\/\//g' > $@

gridscalar.F90: gridscalartemplate.ftn $(LOCALDEPS) $(GLOBALDEPS) $(HEADS)
	sed 's/\([^\#]\)\(\/\/\)/\1\^\^/g' $< | $(FPP) $(FPPFLAGS) -D'_TAG=' \
	-D'__XMAP(_X)=(_X)' -D'__XIMAP(_X)=(_X)' -D'__XGEZERO' -D'__XMAPGEZERO' $(DEFS) $(LOCALDEFS) -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(MSG) $(INCLUDES) -P | sed 's/\^\^/\/\//g' > $@

gridscalar_beta.F90: gridscalartemplate.ftn $(LOCALDEPS) $(GLOBALDEPS) $(HEADS)
	sed 's/\([^\#]\)\(\/\/\)/\1\^\^/g' $< | $(FPP) $(FPPFLAGS) -D'_TAG=_beta' \
	-D'__XMAP(_X)=(one/(kboltz*(_X)))' -D'__XIMAP(_X)=(one/(kboltz*(_X)))' -D'__XGEZERO' -D'__XMAPGEZERO' $(DEFS) $(LOCALDEFS) -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(MSG) $(INCLUDES) -P | sed 's/\^\^/\/\//g' > $@

