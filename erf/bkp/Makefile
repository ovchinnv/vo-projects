DEFS=-D int=integer -D float=real*8 -D bool=logical

FC=ifort
FFLAGS= -O3 -fast -c $(DEFS)

#FC=gfortran
#FFLAGS=-O3 -c # -fbounds-check -frecord-marker=4 -D int=integer -D float=real*8 -D bool=logical $(DEFS)


CC=gcc
CFLAGS=-O3 -c

LD=$(FC)
LDFLAGS=-lm
#LD=$(CC)

all: test.o
	$(LD) $(LDFLAGS) test.o erfsun_2_f90.o
test.o: erfsun_2_f90.o test.F90
	$(FC) $(FFLAGS) test.F90
erf1.o: erf1.c
	$(CC) $(CFLAGS) erf1.c
erfsun_2_f90.o: erfsun_2_f90.c
	$(CC) $(CFLAGS) erfsun_2_f90.c
clean: 
	rm -f *.o