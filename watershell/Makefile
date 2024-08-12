NAME=watershell
ROOT=..
HEADS=watershell.h
OBJS= watershell.o watershell_plugin.o ugrid.o watershell_sbf_spline.o watershell2.o # watershell2.o needs to be last because it redefines _N

DEPS=files dynamol continua output parser constants erf lu unitcell bestfit ppoly

LIB=libwatershell.so # libwatershell.a
LIBS+=$(LIB)

include version.inc
ifeq ($(version),2)
 LOCALMACROS+=-D__USE_VERSION_2
endif

include linkobjs.inc
include $(ROOT)/Makefile.inc
ifeq ($(strip($ACEMD_PLUGIN)),1)
 OBJS+=watershell_plugin_acemd.o
endif

LOCALMACROS+=-include source.defs -D__BINDC -include $(ROOT)/erf/erf_as.def -include $(ROOT)/erf/exp_as.def  #-D__erf=erfo7
LOCALDEPS+=source.defs $(ROOT)/erf/erf_as.def $(ROOT)/erf/exp_as.def

LOCALDEFS+=-fprotect-parens

default : plugin
plugin : $(LIB)
$(LIB) : $(OBJS) $(COBJS) $(LOBJS) linkobjs.inc

$(SOURCES) : $(LOCALDEPS)

watershell.o watershell2.o : ugrid.o

watershell2.o : watershell_sbf_spline.o

watershell_plugin.o : watershell.o watershell2.o

watershell.F90 watershell_plugin.F90 : parser/parser_macro2.src parser/parser_macro3.src watershell_rmft_aux1.ftn

watershell2.F90 : parser/parser_macro2.src parser/parser_macro3.src

watershell_plugin_acemd.o : aceplug.h

aceplug.h:
	wget raw.githubusercontent.com/Acellera/plumed1/master/plumed-1.3/ACEMD/aceplug.h

watershell.F90:watershell.ftn $(LOCALDEPS) $(GLOBALDEPS) $(HEADS)
	@echo " "$<" ==> "$@" ..."
	sed 's/\([^\#]\)\(\/\/\)/\1\^\^/g' $< | $(FPP) $(FPPFLAGS) $(DEFS) $(LOCALDEFS) -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(MSG) $(INCLUDES) -P | sed 's/\^\^/\/\//g'|sed 's/e \?- \?/e-/g'|sed 's/d \?- \?/d-/g'> $@

watershell2.F90:watershell2.ftn $(LOCALDEPS) $(GLOBALDEPS) $(HEADS)
	@echo " "$<" ==> "$@" ..."
	sed 's/\([^\#]\)\(\/\/\)/\1\^\^/g' $< | $(FPP) $(FPPFLAGS) $(DEFS) $(LOCALDEFS) -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(MSG) $(INCLUDES) -P | sed 's/\^\^/\/\//g'|sed 's/e \?- \?/e-/g'|sed 's/d \?- \?/d-/g'> $@

