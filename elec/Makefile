FC=gfortran
FCFLAGS=-c -O3 -fimplicit-none -ffixed-line-length-none -ffree-line-length-none # -fcheck=all -g
FPP=cpp
#FPPFLAGS=-D_FILTER=_POLY3
FPPFLAGS=-D_FILTER=_GAUSS
FPPFLAGS+=-D_NORMALIZE_FILTER

LD=$(FC)
LDFLAGS=

FILES=vars setup run filt3 direct3 direct3_nosym util
OBJS=$(addsuffix .o, $(FILES))
TESTOBJS:=$(filter-out run.o, $(OBJS)) test.o
SRC=$(addsuffix .f90, $(FILES))
DEFS=source.def filters.def
INCLUDES:=$(addprefix -include , $(DEFS))

EXE=elec

default : $(EXE)
$(EXE) : $(OBJS)
	$(LD) $(OBJS) -o $@

$(SRC) : $(DEFS) Makefile

%.o : %.f90
	$(FC) $(FCFLAGS) $< -o $@
%.f90 : %.F90
	$(FPP) $(FPPFLAGS) -P $(INCLUDES) $< > $@
	
clean :
	rm -f $(OBJS) *.mod *.f90 test.o
	
$(OBJS) : $(SRC) vars.o Makefile

test : $(TESTOBJS)
	$(LD) $(TESTOBJS) -o test
	./test
	