ROOT=..
OBJS=erfappx.o erfsun_ref.o erfsun_value.o erf_ext.o
NAME=erf

LOCALMACROS=-D__F2003

include $(ROOT)/Makefile.inc

clean::
	rm -f $(ROOT)/$(NAME).stamp $(SOURCES) $(OBJS)

test.o: erfsun_2_f90.o test.F90
	$(FC) $(FFLAGS) test.F90
