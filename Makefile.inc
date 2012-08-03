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
MPI_HOME=/usr/local/mpich2-intel
#FC=gfortran
#FFLAGS=-O3 -c -frecord-marker=4 -ffixed-line-length-none -ffree-line-length-none -fno-math-errno -fbounds-check -frange-check #-fopenmp 
#FFLAGS=-O3 -c -frecord-marker=4 -ffixed-line-length-none -ffree-line-length-none -fno-math-errno #-fbounds-check -frange-check 
#MPI_HOME=/usr/local/mpich2-gnu
#FC=pathf95
#FFLAGS=-Ofast -fno-math-errno -extend-source -OPT:Olimit=8753 -INLINE -ipa -c
#FFLAGS=-Ofast -fno-math-errno -extend-source -ffortran-bounds-check -OPT:Olimit=8753 -INLINE -ipa -c
#FFLAGS=-extend-source -OPT:Olimit=8753 -INLINE -c
#MPI_HOME=/usr/local/mpich2-path64
#
DEFS=-D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4\
     -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP

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
LDFLAGS=-static -Bstatic -ipo
#LDFLAGS=-static -ipa

DYNAMOL=objects/molecule/dynamol
# root must be defined elswhere
DYNA=$(ROOT)/$(DYNAMOL)

GLOBALMACROS=-include $(ROOT)/source.defs 
MSG=-include $(ROOT)/source.msg
CHMSG=-include $(ROOT)/charmm.msg
TESTMSG=-include $(ROOT)/test.msg
GLOBALDEPS=$(ROOT)/source.defs $(ROOT)/source.msg

MPI_INCLUDE=-I $(MPI_HOME)/include
MPI_LIB=-L $(MPI_HOME)/lib

INCLUDES+=-I $(ROOT)/include
LIBS+=-L $(ROOT)/lib
OBJDIR=$(ROOT)/obj
LOBJS=$(addprefix $(OBJDIR)/,$(LINKOBJS))
DEPS+=$(addprefix $(NAME)/,$(SUBDIRS))

SOURCES=$(OBJS:%.o=%.F90)
CHARMM_SOURCES=$(OBJS:%.o=%.src)

VPATH=$(ROOT) #this is necessary to prevent recompilation under %.stamp rule

FPP=cpp -P
FPPFLAGS = -C -I $(ROOT)

.SUFFIXES:
.SUFFIXES: .o .F90 .a .H .h .c .f .F .stamp

default: $(ROOT)/$(NAME).stamp
all : default
########################### OBJECTS ##########################
%.o:%.f $(LIBDEPS)
	$(FC) $(FFLAGS) $(DEFS) $(LOCALDEFS) $(INCLUDES) $(MPI_INCLUDE) $< -o $@
%.o:%.c $(CLIBDEPS)
	$(CC) $(CFLAGS) $< -o $@
%.o:%.F90 $(LIBDEPS)
	$(FC) $(FFLAGS) $(DEFS) $(LOCALDEFS) $(INCLUDES) $(MPI_INCLUDE) $< -o $@
%.o:%.F $(LIBDEPS)
	$(FC) $(FFLAGS) $(DEFS) $(LOCALDEFS) $(INCLUDES) $(MPI_INCLUDE) $< -o $@
%.a:
	$(AR) $(AFLAGS) $@ $(OBJS) $(LOBJS)
	file=$@ ; if [ -f $$file ] ; then ln -fs `pwd`/$$file ${ROOT}/lib/$$file; 2>/dev/null; fi
###############################################################
# DISTRIBUTION DIRS
DISTRIB=$(addprefix $(ROOT)/, include obj lib)
$(DISTRIB):
	if [ ! -d $@ ] ; then mkdir $@ ; fi

###################### STAMPS KEEP TRACK OF DEPENDENCIES ON DIRECTORIES (OR OTHER COMPOUND OBJECTS): #####
# try searching for module file in uppercase for pathscale
$(ROOT)/$(NAME).stamp:: $(OBJS) $(DISTRIB)
	for file in $(OBJS); do \
	ln -fs `pwd`/$$file ${ROOT}/obj/$$file; \
	file=$${file%.*}.mod ; if [ -f $$file  ] ; then ln -fs `pwd`/$$file ${ROOT}/include/$$file; fi ;\
	file=`echo $${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
	if [ -f $$file  ] ; then ln -fs `pwd`/$$file ${ROOT}/include/$$file; fi ;\
	done; 
	echo -n > $@
####################### GENERIC STAMP RESULTS IN DIFFERENT BEHAVIOR : ####################################
%.stamp:
	DIR=$@ ; DIR=$(ROOT)/$${DIR%.*} ; if [ -d $$DIR ] ; then cd $$DIR; $(MAKE) ; fi

########################## SOURCE CODE #######################
%.h:%.H $(LOCALDEPS) $(GLOBALDEPS)
	$(FPP) $(FPPFLAGS) $(DEFS) $(LOCALDEFS) $< > $@
%.F90:%.ftn $(LOCALDEPS) $(GLOBALDEPS) $(ROOT)/source.msg
	$(FPP) $(FPPFLAGS) $(DEFS) -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(MSG) -P $< > $@
%.src:%.ftn $(LOCALDEPS) $(GLOBALDEPS) $(ROOT)/charmm.msg
	$(FPP) $(FPPFLAGS)         -D__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(CHMSG) -P $< > $@
########################## TEST PROGRAM ######################
test.F90:test.ftn $(LOCALDEPS) $(GLOBALDEPS) $(ROOT)/test.msg
	$(FPP) $(FPPFLAGS) $(DEFS) -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(TESTMSG) -P $< > $@
test :: $(OBJS) test.o force_look
	$(LD) $(LDFLAGS) test.o $(OBJS) $(LOBJS) $(LIBS) -o test
	echo "==================================================="
	echo "Running test :"
	echo "==================================================="
	./test 2>&1 > test.out ; cat test.out
	echo "==================================================="
	echo
######################### OTHER RULES ########################
source: $(SOURCES)
$(OBJS):$(addsuffix .stamp, $(DEPS))

####################### EXECUTABLE ######################
$(EXE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LOBJS) $(LIBS) -o $(EXE)

dynamol:
	$(MAKE) -C $(ROOT)/dynamol
string:
	$(MAKE) -C $(ROOT)/string
chest:
	$(MAKE) ches -C $(ROOT)/chest


force_look:
	true
mod:
	/bin/sh -c 'for file in `ls *.mod *.h 2>/dev/null`; do  ln -fs `pwd`/$$file ${ROOT}/include/$$file; done'
modules: $(OBJS) mod
	for d in $(DIRS); do (cd $$d; $(MAKE) modules); done
linklibs:
	/bin/sh -c 'for file in `ls *.a 2>/dev/null`; do  ln -fs `pwd`/$$file ${ROOT}/lib/$$file; done'
libs : libraries
libraries: $(ROOT)/$(NAME).stamp linklibs
	for d in $(DIRS); do (cd $$d; $(MAKE) libraries); done
.SILENT : %.stamp clean $(DISTRIB) test
clean::
	rm -f *.mod *.a $(ROOT)/$(NAME).stamp $(CHARMM_SOURCES) $(SOURCES) $(OBJS) $(EXE) test.F90 test.dat test.out test.o 
	if [ -f test ] ; then rm -f test ; fi

charmm:
