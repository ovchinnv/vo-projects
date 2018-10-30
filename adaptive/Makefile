NAME=adaptive
ROOT=..
OBJS=adaptive_kforce.o

DEPS=output constants griddata parser files random string/smcv

LOCALMACROS+=-include source.defs
LOCALDEPS+=source.defs

# for linking test executable
LINKOBJS = constants.o rng.o clcg.o output.o parser.o prm.o cmd.o files.o ivector.o isort.o

include $(ROOT)/Makefile.inc
$(SOURCES) : $(LOCALDEPS)

adaptive_kforce.F90: adaptive_template.ftn $(LOCALDEPS) $(GLOBALDEPS) $(HEADS) Makefile
	sed 's/\([^\#]\)\(\/\/\)/\1\^\^/g' $< | $(FPP) $(FPPFLAGS) -D'_NAME=kforce' -D'_MAPTYPE=_log' \
	-D'__XMAP(_X)=log(_X)' -D'__DXMAP(_X)=(one/(_X))' -D'__XIMAP(_X)=exp(_X)' \
	-U__CUSTOM_DIFFUSION -D'__D(_X)=(1d0 + _X)' -D'__DP(_X)=(1d0)' \
	-D'__FORCE(__DUAL, __PAR)=-(__DUAL-__CONCAT3(gridscalar,_MAPTYPE,_lookup_average_value(__PAR)))' \
	-U__EVOLVE_MAPPED \
	-D'__XGTZERO' $(DEFS) $(LOCALDEFS) -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(MSG) $(INCLUDES) -P | sed 's/\^\^/\/\//g' > $@
#	-D'__D(_X)=(1d0)' -D'__DP(_X)=(0d0)' \
#	-D'__D(_X)=(_X)' -D'__DP(_X)=(1d0)' \
