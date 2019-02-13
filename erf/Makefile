ROOT=..
OBJS=erfappx.o erfsun_ref.o erfsun_value.o erf_ext.o
NAME=erf

LOCALMACROS=-D__F2003 -include erf_as.def
LOCALDEPS+=erf_as.def
TESTOBJS=$(OBJS)

include $(ROOT)/Makefile.inc

clean::
	rm -f $(ROOT)/$(NAME).stamp $(SOURCES) $(OBJS)

test.o: test.F90 erf_as.def
	$(FC) $(FFLAGS) test.F90
