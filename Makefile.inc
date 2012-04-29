SHELL=/bin/sh
#
CC=gcc
CFLAGS =-O3 -fno-math-errno -c
#
FC=ifort
FFLAGS=-O3 -fast -finline -fp-model fast=2 -xhost -ipo -fno-math-errno -no-inline-factor -no-inline-min-size -c
#FFLAGS= -O3 -c -check-all -extend-source -fast -fno-math-errno -finline -fp-model fast=2 -xhost -ipo  #-ipo 
#FFLAGS= -O3 -c -extend-source -fast -fno-math-errno -finline -fp-model fast=2 -xSSE3 -ipo # -no-inline-factor -no-inline-min-size #-ipo 
#FFLAGS= -O3 -c -extend-source -fast -fno-math-errno -finline -fp-model fast=2 -xhost -no-ipo # -no-inline-factor -no-inline-min-size #-ipo 
#FC=gfortran
#FFLAGS=-O3 -c -frecord-marker=4 -ffixed-line-length-132 -ffree-line-length-132 -fno-math-errno -fbounds-check -frange-check #-fopenmp 
#FFLAGS=-O3 -c -frecord-marker=4 -ffixed-line-length-132 -ffree-line-length-132 -fno-math-errno #-fbounds-check -frange-check 
#
DEFS=-D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4
# ! note: the VEC routines are probably slower

# note : for single precision define
#
#DEFS=$(DEFS) -D DEBUG
#
# preprocess with C++ style preprocessor
# -C flag preserves comments

AR=xiar
#AR=ar
AFLAGS=rv

LD=$(FC)
LDFLAGS=-static -Bstatic

DYNAMOL=objects/molecule/dynamol
# root must be defined elswhere
DYNA=$(ROOT)/$(DYNAMOL)

GLOBALMACROS=-include $(ROOT)/source.defs 
MSG=-include $(ROOT)/source.msg
CHMSG=-include $(ROOT)/charmm.msg
GLOBALDEPS=$(ROOT)/source.defs $(ROOT)/source.msg

MPI_HOME=/usr/local/mpich2-intel

MPI_INCLUDE=-I $(MPI_HOME)/include
MPI_LIB=-L $(MPI_HOME)/lib

SOURCES=$(OBJS:%.o=%.F90)
CHARMM_SOURCES=$(OBJS:%.o=%.src)

FPP=cpp -P
FPPFLAGS = -C

.SUFFIXES:
.SUFFIXES: .o .F90 .a .H .h .c .f .F

default : all
all: $(OBJS)
########################### OBJECTS ##########################
%.o:%.F90 $(LIBDEPS)
	$(FC) $(FFLAGS) $(DEFS) $(LOCALDEFS) $(INCLUDES) $(MPI_INCLUDE) $<
%.o:%.f $(LIBDEPS)
	$(FC) $(FFLAGS) $(DEFS) $(LOCALDEFS) $(INCLUDES) $(MPI_INCLUDE) $<
%.o:%.F $(LIBDEPS)
	$(FC) $(FFLAGS) $(DEFS) $(LOCALDEFS) $(INCLUDES) $(MPI_INCLUDE) $<
%.o:%.c $(CLIBDEPS)
	$(CC) $(CFLAGS) $<
%.a:
	$(AR) $(AFLAGS) $@ $(OBJS)
########################## SOURCE CODE #######################
%.h:%.H $(LOCALDEPS) $(GLOBALDEPS)
	$(FPP) $(FPPFLAGS) $(DEFS) $(LOCALDEFS) $< > $@
%.F90:%.ftn $(LOCALDEPS) $(GLOBALDEPS)
	$(FPP) $(FPPFLAGS) $(DEFS) -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(MSG) -P $< > $@
%.src:%.ftn $(LOCALDEPS) $(GLOBALDEPS) 
	$(FPP) $(FPPFLAGS)         -D__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(CHMSG) -P $< > $@
######################### OTHER RULES ########################
source: $(SOURCES)
#$(OBJS):$(SOURCES)
$(DIRS): force_look
	$(MAKE) -C $@ modules

output:
	$(MAKE) -C $(ROOT)/output modules
parser:
	$(MAKE) -C $(ROOT)/parser modules
ivector:
	$(MAKE) -C $(ROOT)/vectors modules
vectors:
	$(MAKE) -C $(ROOT)/vectors modules
files:
	$(MAKE) -C $(ROOT)/files modules
constants:
	$(MAKE) -C $(ROOT)/constants modules


force_look:
	true
mod:
	/bin/sh -c 'for file in `ls *.mod 2>/dev/null`; do  ln -fs `pwd`/$$file ${ROOT}/modules/$$file; done'
modules: $(OBJS) mod
	for d in $(DIRS); do (cd $$d; $(MAKE) modules); done
lib: force_look
	/bin/sh -c 'for file in `ls *.a 2>/dev/null`; do  ln -fs `pwd`/$$file ${ROOT}/lib/$$file; done'
libraries: $(OBJS) lib
	for d in $(DIRS); do (cd $$d; $(MAKE) libraries); done
.SILENT : force_look clean modules mod
clean: 
	rm -f *.mod $(CHARMM_SOURCES) $(SOURCES) $(OBJS) $(EXE)

.PHONY : parser output vectors files constants clean
