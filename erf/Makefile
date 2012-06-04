ROOT=..
OBJS=erfappx.o
NAME=erf


include $(ROOT)/Makefile.inc

$(ROOT)/$(NAME).stamp :: erfsun_2_f90.o
erfsun_2_f90.o : erfsun_2_f90.c

clean::
	rm -f $(ROOT)/$(NAME).stamp $(SOURCES) $(OBJS) erfsun_2_f90.o
