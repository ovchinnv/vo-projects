ROOT=..
OBJS=enm.o
NAME=enm
DEPS=output lu constants files multicom parselist parser

LINKOBJS=$(addsuffix .o, $(DEPS)) ivector.o multicom_aux.o

include $(ROOT)/Makefile.inc

fit :: $(OBJS) fit.o
	$(LD) $(LDFLAGS) fit.o $(OBJS) $(LOBJS) $(LIBS) $(MPI_LIB) -o fit

fit.o: $(OBJS) fit.F90
fit.F90 : fit.ftn

clean :: 
	rm -f fit fit.o fit.F90
