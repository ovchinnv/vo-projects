ROOT=..
OBJS=erfappx.o erfsun_ref.o erfsun_value.o erf_ext.o ierfappx.o
HEADS=erf.h
NAME=erf

LOCALMACROS=-D__F2003 -include erf_as.def -include exp_as.def
LOCALDEPS+=erf_as.def exp_as.def
FPPFLAGS+=
TESTOBJS=$(OBJS)

include $(ROOT)/Makefile.inc
LOCALDEFS+=-fprotect-parens

clean::
	rm -f $(ROOT)/$(NAME).stamp $(SOURCES) $(OBJS)

test.F90:test.ftn $(LOCALDEPS) $(GLOBALDEPS) $(ROOT)/test.msg
	$(FPP) $(FPPFLAGS) $(DEFS) $(LOCALDEFS) -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(TESTMSG) $(DEBUG) -P $< |sed 's/e \?- \?/e-/g'|sed 's/d \?- \?/d-/g'>$@

test.o: test.F90 erf_as.def
	$(FC) $(FFLAGS) test.F90

ierfappx.F90:ierfappx.ftn $(LOCALDEPS) $(GLOBALDEPS) $(HEADS)
	@echo " "$<" ==> "$@" ..."
	$(FPP) $(FPPFLAGS) $(DEFS) $(LOCALDEFS) -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(MSG) $(INCLUDES) -P $< |sed 's/e \?- \?/e-/g'|sed 's/d \?- \?/d-/g'>$@
