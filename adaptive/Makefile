NAME=adaptive
ROOT=..
OBJS=adaptive_kforce.o

DEPS=output constants griddata parser files random

LOCALMACROS+=-include source.defs
LOCALDEPS+=source.defs

include $(ROOT)/Makefile.inc
$(SOURCES) : $(LOCALDEPS)

adaptive_kforce.F90: adaptive_template.ftn $(LOCALDEPS) $(GLOBALDEPS) $(HEADS)
	sed 's/\([^\#]\)\(\/\/\)/\1\^\^/g' $< | $(FPP) $(FPPFLAGS) -D'_NAME=kforce' -D'_MAPTYPE=_log' \
	-D'__XMAP(_X)=log(_X)' -D'__DXMAP(_X)=(one/(_X))' -D'__XIMAP(_X)=exp(_X)' -D'__DIFF(_X)=1' \
	-D'__DRIFT(_DUAL,_PAR)=-bet*(_DUAL-__CONCAT3(gridscalar,_MAPTYPE,_lookup_average_value(_PAR)))' \
	-D'__XGTZERO' $(DEFS) $(LOCALDEFS) -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(MSG) $(INCLUDES) -P | sed 's/\^\^/\/\//g' > $@
