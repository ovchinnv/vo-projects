FC=gfortran
FCFLAGS=-c -O3 -fimplicit-none -ffixed-line-length-none -ffree-line-length-none -g
FPP=cpp
FPPFLAGS=-D_FILTER=_POLY3
INCLUDES=-include source.def
LD=$(FC)
LDFLAGS=

FILES=vars setup run filt3p direct3p direct3p_nosym util
OBJS=$(addsuffix .o, $(FILES))
SRC=$(addsuffix .f90, $(FILES))
EXE=elec

default : $(EXE)
$(EXE) : $(OBJS)
	$(LD) $(OBJS) -o $@

$(SRC) : source.def

%.o : %.f90
	$(FC) $(FCFLAGS) $< -o $@
%.f90 : %.F90
	$(FPP) $(FPPFLAGS) -P $(INCLUDES) $< > $@
	
clean :
	rm -f $(OBJS) *.mod *.f90
	
$(OBJS) : vars.o Makefile