SHELL=/bin/sh
#
CC=gcc
CFLAGS =-O3 -fno-math-errno -c -fPIC
#---INTEL
#FC=ifort
#FFLAGS=-O2 -c
#FFLAGS=-O3 -fast -finline -fp-model fast=2 -xhost -ipo -fno-math-errno -no-inline-factor -no-inline-min-size -c -g
#FFLAGS= -O3 -fast -c -ipo -debug#-fast -fno-math-errno -finline -fp-model fast=2 -xhost -no-ipo #-ipo 
#FFLAGS= -O0 -c -g -check all -debug extended -debug-parameters all -traceback -no-ipo #-fast -fno-math-errno -finline -fp-model fast=2 -xhost -no-ipo #-ipo 
#FFLAGS= -O3 -c -fast -fno-math-errno -finline -fp-model fast=2 -xSSE3 -ipo # -no-inline-factor -no-inline-min-size #-ipo #xSSE3 for pentium
#FFLAGS= -O2 -c -fast -fno-math-errno -finline -fp-model fast=2 -xhost -no-ipo # -no-inline-factor -no-inline-min-size #-ipo 
#MPI_HOME=/usr/local/mpich2-intel
#---GNU
FC=gfortran
#FFLAGS=-O3 -c -g -frecord-marker=4 -ffixed-line-length-none -ffree-line-length-none -fno-math-errno -Wall  #-fopenmp 
#FFLAGS=-O3 -c -g -frecord-marker=4 -ffixed-line-length-none -ffree-line-length-none -fno-math-errno -fcheck=all  -Wall  #-fopenmp 
#FFLAGS=-O3 -c -frecord-marker=4 -ffixed-line-length-none -ffree-line-length-none -fcheck=all
FFLAGS=-O3 -c -frecord-marker=4 -ffixed-line-length-none -ffree-line-length-none -fno-math-errno
#MPI_HOME=/usr/local/mpich2-gnu
MPI_HOME=/usr/lib/openmpi
OMP=-D__OMP -fopenmp
LOMP=-lgomp
#---PATHSCALE
#FC=pathf95
#FFLAGS=-Ofast -fno-math-errno -extend-source -INLINE -c -ipa
#FFLAGS=-Ofast -fno-math-errno -extend-source -ffortran-bounds-check -OPT:Olimit=8753 -INLINE -ipa -c
#FFLAGS=-extend-source -OPT:Olimit=8753 -INLINE -c
#MPI_HOME=/usr/local/mpich2-path64
#FC=/usr/bin/g95
#FFLAGS=-O3 -c -g -ffree-line-length-huge -ffixed-line-length-132 -fno-second-underscore -cpp
#MPI_HOME=/usr/local/mpich2-g95
#
# G95 reguires no space between "-D" and the name
FORTRAN_COMPILER=$(notdir $(FC))
DEFS+=-Dint=integer -Dfloat=real*8 -Dbool=logical -Dmpiint=MPI_INTEGER -Dmpifloat=MPI_REAL8 -Dmpichar=MPI_CHARACTER -Dmpiint4=MPI_INTEGER4\
     -Dmpibool=MPI_LOGICAL -Dint1="integer(kind=KIND('a'))" -D__RCOMP -D__DMOL -D__CHES -D__EDDIE -D__$(FORTRAN_COMPILER) \
     -Dint4mpi=integer -D__TIMER  -D__PARALLEL #-D__DEBUG -D__PATHSCALE #-Dint4mpi="integer(KIND=MPI_INTEGER_KIND)" 

STRING_FLAGS=-D__STRING_SERIAL -D__ACEMD -D__NAMD $(OMP)
DEFS+=$(STRING_FLAGS)
#
# preprocess with C++ style preprocessor
# -C flag preserves comments

AR=ar
AFLAGS=rvcs

LD=$(FC)
#LDFLAGS=-no-ipo -static -Bstatic #-ipo
#LDFLAGS=-static -Bstatic -ipo
LDFLAGS=-dynamic -Bdynamic $(LOMP) #-Wl,-z,defs # to list unresolved symbols during linking
#LDFLAGS=-static -ipa

#EXTRA_LIBS+=-lgfortran
# extra libs for ifort
#EXTRA_LIBS+=-lirc -limf -lifcore -ldl
#
# G95
#EXTRA_LIBS+=-lf95 -lgcc

EXTRA_LIBS+=-llapack

PIC=-fPIC
GLOBALMACROS=-include $(ROOT)/source.defs
GLOBALMACROS+=-include $(ROOT)/contexts.def

MSG=-include $(ROOT)/source.msg
CHMSG=-include $(ROOT)/charmm.msg
TESTMSG=-include $(ROOT)/test.msg
GLOBALDEPS=$(ROOT)/source.defs $(ROOT)/source.msg 
GLOBALDEPS+=$(ROOT)/contexts.def $(ROOT)/setcomm.def  $(ROOT)/exec.def

MPI_INCLUDE=-I $(MPI_HOME) # -I $(MPI_HOME)/include -I $(MPI_HOME)/lib
MPI_LIBS=-L $(MPI_HOME) -L $(MPI_HOME)/lib -lpthread -lmpich -lmpichf90 -lmpl # -lmpi -lmpi_f90 -lmpi_f77 # mpich
MPI_LIBS=-L $(MPI_HOME) -L $(MPI_HOME)/lib -lpthread -lmpi -lmpi_f90 -lmpi_f77 # older openmpi
MPI_LIBS=-L $(MPI_HOME) -lmpi -lmpi_mpifh -lmpi_usempif08  # newer openmpi w/ F08
#MPI_INCLUDE=
#MPI_LIBS=

INCLUDES+=-I $(ROOT)/include
LIBS+=-L $(ROOT)/lib
LIBS+=$(EXTRA_LIBS)
OBJDIR=$(ROOT)/obj
LOBJS=$(addprefix $(OBJDIR)/,$(LINKOBJS))
DEPS+=$(addprefix $(NAME)/,$(SUBDIRS))

SOURCES=$(OBJS:%.o=%.F90)
CHARMM_SOURCES=$(OBJS:%.o=%.src)
CHARMM_MAJOR_VERSION=37

VPATH=$(ROOT) #this is necessary to prevent recompilation under %.stamp rule

FPP=cpp
FPPFLAGS = -P -C -I $(ROOT) -nostdinc # -nostdinc : avoid including files automatically

DEBUG = -g

.SUFFIXES:
.SUFFIXES: .o .F90 .a .so .H .h .c .f .F .stamp

all : default
default: $(ROOT)/$(NAME).stamp
########################### OBJECTS ##########################
%.o:%.f $(LIBDEPS)
	$(FC) $(FFLAGS) $(DEFS) $(PIC) $(LOCALDEFS) $(INCLUDES) $(MPI_INCLUDE) $< -o $@
%.o:%.c $(CLIBDEPS)
	$(CC) $(LOCALDEFS) $(CFLAGS) $(PIC) $< -o $@
%.o:%.F90 $(LIBDEPS)
	$(FC) $(FFLAGS) $(DEFS) $(PIC) $(LOCALDEFS) $(INCLUDES) $(MPI_INCLUDE) $< -o $@
%.mod:%.F90 $(LIBDEPS)
	$(FC) $(FFLAGS) $(DEFS) $(PIC) $(LOCALDEFS) $(INCLUDES) $(MPI_INCLUDE) $< -o /dev/null
	file=$@ ; if [ -f $$file ] ; then ln -fs `pwd`/$$file ${ROOT}/include/$$file; 2>/dev/null; fi
%__.o:%_.F90 $(LIBDEPS)
	$(FC) $(FFLAGS) $(DEFS) $(PIC) $(LOCALDEFS) $(INCLUDES) $(MPI_INCLUDE) $< -o $@
%.o:%.F $(LIBDEPS)
	$(FC) $(FFLAGS) $(DEFS) $(PIC) $(LOCALDEFS) $(INCLUDES) $(MPI_INCLUDE) $< -o $@
%.a:    $(OBJS) $(HEADS)
	$(AR) $(AFLAGS) $@ $(LOBJS)
	file=$@ ; if [ -f $$file ] ; then ln -fs `pwd`/$$file ${ROOT}/lib/$$file; 2>/dev/null; fi
%.so:   $(OBJS) $(HEADS)
	$(LD) $(LDFLAGS) -shared $(LOBJS) -o $@
	file=$@ ; if [ -f $$file ] ; then ln -fs `pwd`/$$file ${ROOT}/lib/$$file; 2>/dev/null; fi
###############################################################
# DISTRIBUTION DIRS
DISTRIB=$(addprefix $(ROOT)/, include obj lib)
$(DISTRIB):
	if [ ! -d "$@" ] ; then mkdir -p $@ ; fi ; #echo $@

###################### STAMPS KEEP TRACK OF DEPENDENCIES ON DIRECTORIES (OR OTHER COMPOUND OBJECTS): #####
# try searching for module file in uppercase for pathscale
$(ROOT)/$(NAME).stamp:: $(OBJS) $(MODS) $(HEADS) $(DISTRIB)
	for file in $(OBJS); do \
	ln -fs `pwd`/$$file ${ROOT}/obj/$$file; \
	file=$${file%.*}.mod ; if [ -f $$file  ] ; then ln -fs `pwd`/$$file ${ROOT}/include/$$file; fi ;\
	file=`echo $${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
	if [ -f $$file  ] ; then ln -fs `pwd`/$$file ${ROOT}/include/$$file; fi ;\
	done; 
	for file in $(HEADS); do \
	ln -fs `pwd`/$$file ${ROOT}/include/$$file; \
        done;
	touch $@
####################### GENERIC STAMP RESULTS IN DIFFERENT BEHAVIOR : ####################################
%.stamp::
	DIR=$@ ; DIR=$(ROOT)/$${DIR%.*} ; echo $$DIR ; if [ -d $$DIR ] ; then cd $$DIR; $(MAKE) ; fi
########################## SOURCE CODE #######################
%.h:%.H $(LOCALDEPS) $(GLOBALDEPS)
	$(FPP) $(FPPFLAGS) $(DEFS) $(LOCALDEFS) $(GLOBALMACROS) $(LOCALMACROS) -P $< > $@
%.F90:%.ftn $(LOCALDEPS) $(GLOBALDEPS)
	sed 's/\([^\#]\)\(\/\/\)/\1\^\^/g' $< | $(FPP) $(FPPFLAGS) $(DEFS) $(LOCALDEFS) -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(MSG) $(INCLUDES) -P | sed 's/\^\^/\/\//g' > $@
%.src:%.ftn $(LOCALDEPS) $(GLOBALDEPS) $(ROOT)/charmm.msg
	sed 's/\([^\#]\)\(\/\/\)/\1\^\^/g' $< | $(FPP) $(FPPFLAGS) $(STRING_FLAGS) -D__CHARMM -D__CHARMM_VER=$(CHARMM_MAJOR_VERSION) $(GLOBALMACROS) $(LOCALMACROS) $(CHMSG) -P | sed 's/\^\^/\/\//g' > $@
%.ch:%.H $(LOCALDEPS) $(GLOBALDEPS)
	$(FPP) $(FPPFLAGS)  $(STRING_FLAGS) -D__CHARMM -D__CHARMM_VER=$(CHARMM_MAJOR_VERSION) $(CHMSG) $(GLOBALMACROS) $(LOCALMACROS) -P $< > $@
########################## TEST PROGRAM ######################
test.F90:test.ftn $(LOCALDEPS) $(GLOBALDEPS) $(ROOT)/test.msg
	$(FPP) $(FPPFLAGS) $(DEFS) -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(TESTMSG) $(DEBUG) -P $< > $@
test :: $(OBJS) test.o force_look
	$(LD) $(LDFLAGS) test.o $(OBJS) $(LOBJS) $(LIBS) $(MPI_LIBS) $(DEBUG) -o test
	echo "==================================================="
	echo "Running test :"
	echo "==================================================="
	time ./test > test.out 2>&1 ; cat test.out
	echo "==================================================="
	echo
######################### OTHER RULES ########################
source: $(SOURCES)
$(OBJS):$(addsuffix .stamp, $(DEPS))

####################### EXECUTABLE ######################
$(EXE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LOBJS) $(LIBS) $(MPI_LIBS) -o $(EXE)

dynamol:
	$(MAKE) -C $(ROOT)/dynamol
string:
	$(MAKE) -C $(ROOT)/string
continua:
	$(MAKE) -C $(ROOT)/continua

force_look:
	true
mod:
	/bin/sh -c 'for file in `ls *.mod *.h 2>/dev/null`; do  ln -fs `pwd`/$$file ${ROOT}/include/$$file; done'
modules: $(OBJS) mod
	for d in $(DIRS); do (cd $$d; $(MAKE) modules); done
# LIBRARIES
linklibs: libs mod
	/bin/sh -c 'for file in `ls *.a *.so 2>/dev/null`; do  ln -fs `pwd`/$$file ${ROOT}/lib/$$file; done'
libs : libraries
libraries: $(ROOT)/$(NAME).stamp
	for d in $(DIRS); do (cd $$d; $(MAKE) libraries); done

.SILENT : %.stamp $(DISTRIB) force_look
#.SILENT : clean
clean::
	rm -f *.mod *.a *.so $(ROOT)/$(NAME).stamp $(CHARMM_SOURCES) $(SOURCES) $(OBJS) $(EXE) $(HEADS) test.F90 test.dat test.out test.o 
	if [ -f test ] ; then rm -f test ; fi
	pushd ../ >/dev/null; rm -f $(ROOT)/$(NAME).stamp ; popd >/dev/null
charmm::$(addsuffix .stamp, $(CHMDEPS))
