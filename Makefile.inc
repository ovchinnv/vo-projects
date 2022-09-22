SHELL=/bin/sh
include $(ROOT)/configure.inc
ifndef DEBUG
 DEBUG=0
endif
#
ifndef CC
CC=gcc
endif
ifeq ($(strip $(DEBUG)),1)
CFLAGS =-O0 -g -fno-math-errno -c -fPIC
else
CFLAGS =-O3 -fno-math-errno -c -fPIC
endif
#
ifndef CUDA
 CUDA=0
endif
#
#---GNU
ifndef FC
FC=gfortran
endif

# below, we also treat gfortran-* as gfortran : replace - with space, and take the first word
ifeq ($(word 1, $(subst -, ,$(strip $(FC)))),gfortran)
ifeq ($(strip $(DEBUG)),1)
# WARNING: WHEN DEBUGGING DO NOT USE THE LINE BELOW. OPTIMIZATION DESTROYS SYNCHRONIZATION
#FFLAGS=-O3 -c -frecord-marker=4 -ffixed-line-length-none -ffree-line-length-none -fcheck=all -g
 FFLAGS=-O1 -c -g -frecord-marker=4 -ffixed-line-length-none -ffree-line-length-none -fcheck=all
# DEFS+=-D__DEBUG
else
# NOTE : Ofast math has -fno-protect-parens which causes __erfo7 to fail in some cases
# FFLAGS=-O3 -c -frecord-marker=4 -ffixed-line-length-none -ffree-line-length-none -fno-math-errno #-fno-protect-parens #-Wall -fopenmp
 FFLAGS=-Ofast -c -frecord-marker=4 -ffixed-line-length-none -ffree-line-length-none -fno-math-errno -fno-stack-arrays
endif
#---INTEL
else ifeq ($(strip $(FC)),ifort)
ifeq ($(strip $(DEBUG)),1)
 FFLAGS= -O0 -c -g -check all -debug extended -debug-parameters all -traceback -no-ipo #-fast -fno-math-errno -finline -fp-model fast=2 -xhost -no-ipo #-ipo 
else
#FFLAGS= -O3 -fast -finline -fp-model fast=2 -xhost -ipo -fno-math-errno -no-inline-factor -no-inline-min-size -c -g
#FFLAGS= -O3 -fast -c -ipo -debug#-fast -fno-math-errno -finline -fp-model fast=2 -xhost -no-ipo #-ipo 
#FFLAGS= -O3 -c -fast -fno-math-errno -finline -fp-model fast=2 -xSSE3 -ipo # -no-inline-factor -no-inline-min-size #-ipo #xSSE3 for pentium
 FFLAGS= -O2 -c -fast -fno-math-errno -finline -fp-model fast=2 -xhost -no-ipo # -no-inline-factor -no-inline-min-size #-ipo 
endif
endif

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
# NOTE : G95 reguires no space between "-D" and the name
FORTRAN_COMPILER=$(notdir $(FC))
ifeq ($(strip $(FORTRAN_COMPILER)),pathf90)
DEFS+=_D__PATHSCALE
else ifeq ($(strip $(FORTRAN_COMPILER)),pathf95)
DEFS+=_D__PATHSCALE
else ifeq ($(strip $(FORTRAN_COMPILER)),pathcc)
DEFS+=_D__PATHSCALE
else ifeq ($(strip $(FORTRAN_COMPILER)),pathCC)
DEFS+=_D__PATHSCALE
endif

DEFS+=-Dint=integer -Dbool=logical -Dint1=integer\(kind=KIND\(\'a\'\)\) -D__RCOMP -D__DMOL -D__CHES -D__EDDIE -D__$(FORTRAN_COMPILER) -D__TIMER

ifndef DOUBLE_PRECISION
 DOUBLE_PRECISION=1
endif

ifeq ($(strip $(DOUBLE_PRECISION)),1)
 DEFS+=-Dfloat=real*8
 DEFS+=-Dmpifloat=MPI_REAL8
else
 DEFS+=-Dfloat=real*4
 DEFS+=-Dmpifloat=MPI_REAL4
 LOCALDEFS+=-D__SINGLE
endif

ifndef USE_PARALLEL
 USE_PARALLEL=0
endif

ifeq ($(strip $(USE_PARALLEL)),1)
 DEFS+=-D__PARALLEL \
       -Dmpiint=MPI_INTEGER -Dmpichar=MPI_CHARACTER -Dmpiint4=MPI_INTEGER4 -Dmpiint8=MPI_INTEGER8 -Dmpibool=MPI_LOGICAL \
       -Dint4mpi="integer(KIND=MPI_INTEGER_KIND)"

ifndef MPI_HOME
 MPI_HOME=/usr/lib/openmpi
endif
 MPI_INCLUDE=-I $(MPI_HOME) -I $(MPI_HOME)/include -I $(MPI_HOME)/lib
 MPI_LIBS=-L $(MPI_HOME) -L $(MPI_HOME)/lib -lpthread -lmpich -lmpichf90 -lmpl # -lmpi -lmpi_f90 -lmpi_f77 # mpich
 MPI_LIBS=-L $(MPI_HOME) -L $(MPI_HOME)/lib -lpthread -lmpi -lmpi_f90 -lmpi_f77 # older openmpi
 MPI_LIBS=-L $(MPI_HOME) -lmpi -lmpi_mpifh -lmpi_usempif08  # newer openmpi w/ F08

else
 DEFS+=-Dint4mpi=integer \
       -Dmpiint=MPI_INTEGER -Dmpichar=MPI_CHARACTER -Dmpiint4=MPI_INTEGER4 -Dmpiint8=MPI_INTEGER8 -Dmpibool=MPI_LOGICAL # for now, needed for FTSM compatibility

 MPI_INCLUDE=
 MPI_LIBS=
endif

STRING_FLAGS= -D__ACEMD -D__NAMD $(OMP)
ifndef COMPILE_FOR_CHARMM
COMPILE_FOR_CHARMM=0
endif
ifeq ($(strip $(COMPILE_FOR_CHARMM)),0)
DEFS+=$(STRING_FLAGS)
endif
#
DEFS+=$(EXTRA_DEFINES)
FFLAGS+=$(EXTRA_FC_FLAGS)

# preprocess with C++ style preprocessor
# -C flag preserves comments

AR=ar
AFLAGS=rvcs

LD=$(FC)
#LDFLAGS=-no-ipo -static -Bstatic #-ipo
#LDFLAGS=-static -Bstatic -ipo
#LDFLAGS=-Ofast -dynamic -Bdynamic $(LOMP) -Wl,-z,defs # to list unresolved symbols during linking
LDFLAGS=-Ofast -dynamic -Bdynamic $(LOMP) -Wl,-z,defs # to list unresolved symbols during linking
#LDFLAGS=-static -ipa

ifeq ($(strip $(CUDA)),1)
__CUDA=1
# NVIDIA CUDA COMPILER ######
NVCC=nvcc
#NVCCFLAGS=-O3 -m32 -Xcompiler -mfloat-abi=hard -gencode arch=compute_32,code=sm_32 -ccbin $(CC) --use_fast_math -c
NVCCFLAGS=-O3 -ccbin $(CC) --use_fast_math -c
#-gencode arch=compute_30,code=sm_30 -c
#NVCCFLAGS=-O3 -m32 -Xcompiler -mfloat-abi=hard -gencode arch=compute_32,code=sm_32 -ccbin $(CC) -c
#NVCCFLAGS=-O3 -m32 -Xcompiler -mfloat-abi=hard -gencode arch=compute_32,code=sm_32 -ccbin $(CC) --prec-div=false --ftz=true --prec-sqrt=false -c
CUINCLUDES+=-I/usr/local/cuda/samples/common/inc -I/opt/cuda/samples/common/inc 
EXTRA_LIBS+=-L/usr/local/cuda/lib -L/opt/cuda/lib64  -lcudart -lcuda 
endif
#
EXTRA_LIBS+=-lgfortran -lstdc++ -lgcc -lm
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

# NO MODIFICATIONS SHOULD BE NECESSARY BELOW THIS LINE
#####################################################################################################################
#
INCLUDES+=-I $(ROOT)/include
LIBS+=-L $(ROOT)/lib
LIBS+=$(EXTRA_LIBS)
OBJDIR=$(ROOT)/obj
LOBJS=$(addprefix $(OBJDIR)/,$(LINKOBJS))
DEPS+=$(addprefix $(NAME)/,$(SUBDIRS))

SOURCES=$(OBJS:%.o=%.F90)
CHARMM_SOURCES=$(OBJS:%.o=%.src)
ifndef CHARMM_MAJOR_VERSION
CHARMM_MAJOR_VERSION=37
endif

VPATH=$(ROOT) #this is necessary to prevent recompilation under %.stamp rule

FPP=cpp
FPPFLAGS += -P -C -I $(ROOT) -nostdinc # -nostdinc : avoid including files automatically

DEBUG = -g

.SUFFIXES:
.SUFFIXES: .o .F90 .a .so .H .h .c .cu .f .F .stamp

all : default
default: $(ROOT)/$(NAME).stamp
########################### OBJECTS ##########################
%.o:%.f $(LIBDEPS)
	@echo " "$<" ==> "$@" ..."
	$(FC) $(FFLAGS) $(DEFS) $(PIC) $(LOCALDEFS) $(INCLUDES) $(MPI_INCLUDE) $< -o $@
%.o:%.c $(CLIBDEPS) $(HEADS)
	@echo " "$<" ==> "$@" ..."
	$(CC) $(LOCALDEFS) $(CFLAGS) $(PIC) $(GLOBALMACROS) $(LOCALMACROS) $< -o $@
%.o:%.cu $(HEADS)
	$(NVCC) $(NVCCFLAGS) $(LOCALDEFS) $(CUINCLUDES) $(LOCALMACROS) $< -o $@
%.o:%.F $(LIBDEPS)
	@echo " "$<" ==> "$@" ..."
	$(FC) $(FFLAGS) $(DEFS) $(PIC) $(LOCALDEFS) $(INCLUDES) $(MPI_INCLUDE) $< -o $@
%.o:%.F90 $(LIBDEPS)
	@echo " "$<" ==> "$@" ..."
	$(FC) $(FFLAGS) $(DEFS) $(PIC) $(LOCALDEFS) $(INCLUDES) $(MPI_INCLUDE) $< -o $@
%.mod:%.F90 $(LIBDEPS)
	@echo " "$<" ==> "$@" ..."
	$(FC) $(FFLAGS) $(DEFS) $(PIC) $(LOCALDEFS) $(INCLUDES) $(MPI_INCLUDE) $< -o /dev/null
	file=$@ ; if [ -f $$file ] ; then ln -fs `pwd`/$$file ${ROOT}/include/$$file; 2>/dev/null; fi
%__.o:%_.F90 $(LIBDEPS)
	$(FC) $(FFLAGS) $(DEFS) $(PIC) $(LOCALDEFS) $(INCLUDES) $(MPI_INCLUDE) $< -o $@
%.a:    $(OBJS) $(HEADS)
	$(AR) $(AFLAGS) $@ $(LOBJS) $^
	file=$@ ; if [ -f $$file ] ; then ln -fs `pwd`/$$file ${ROOT}/lib/$$file; 2>/dev/null; fi
%.so:   $(OBJS) $(HEADS)
	$(LD) $(LDFLAGS) -shared $(sort $(LOBJS) $(OBJS) $(COBJS)) -o $@
	file=$@ ; if [ -f $$file ] ; then ln -fs `pwd`/$$file ${ROOT}/lib/$$file; 2>/dev/null; fi
###############################################################
# DISTRIBUTION DIRS
DISTRIB=$(addprefix $(ROOT)/, include obj lib)
$(DISTRIB):
	if [ ! -d "$@" ] ; then mkdir -p $@ ; fi ; #echo $@

###################### STAMPS KEEP TRACK OF DEPENDENCIES ON DIRECTORIES (OR OTHER COMPOUND OBJECTS): #####
# try searching for module file in uppercase for pathscale
$(ROOT)/$(NAME).stamp:: $(OBJS) $(MODS) $(HEADS) $(DISTRIB) $(addprefix $(ROOT)/, $(addsuffix .stamp, $(DEPS)))
	@for file in $(OBJS); do \
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
	@DIR=$@ ; DIR=$${DIR%.*} ; echo $$DIR ; if [ -d $$DIR ] ; then $(MAKE) -C $$DIR ; fi
########################## SOURCE CODE #######################
%.h:%.H $(LOCALDEPS) $(GLOBALDEPS)
	$(FPP) $(FPPFLAGS) $(DEFS) $(LOCALDEFS) $(GLOBALMACROS) $(LOCALMACROS) -P $< > $@
%.hpp:%.H $(LOCALDEPS) $(GLOBALDEPS)
	$(FPP) -D__CPLUSPLUS $(FPPFLAGS) $(DEFS) $(LOCALDEFS) $(GLOBALMACROS) $(LOCALMACROS) -P $< > $@
%.F90:%.ftn $(LOCALDEPS) $(GLOBALDEPS) $(HEADS)
	@echo " "$<" ==> "$@" ..."
	sed 's/\([^\#]\)\(\/\/\)/\1\^\^/g' $< | $(FPP) $(FPPFLAGS) $(DEFS) $(LOCALDEFS) -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(MSG) $(INCLUDES) -P | sed 's/\^\^/\/\//g' > $@
%_.F90:%.ftn $(LOCALDEPS) $(GLOBALDEPS) $(HEADS)
	sed 's/\([^\#]\)\(\/\/\)/\1\^\^/g' $< | $(FPP) $(FPPFLAGS) $(DEFS) $(LOCALDEFS) -D__BINDC -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(MSG) $(INCLUDES) -P | sed 's/\^\^/\/\//g' > $@
%.src:%.ftn $(LOCALDEPS) $(GLOBALDEPS) $(ROOT)/charmm.msg
	sed 's/\([^\#]\)\(\/\/\)/\1\^\^/g' $< | $(FPP) $(FPPFLAGS) -D__CHARMM -D__CHARMM_VER=$(CHARMM_MAJOR_VERSION) $(GLOBALMACROS) $(LOCALMACROS) $(CHMSG) -P | sed 's/\^\^/\/\//g' > $@
%.ch:%.H $(LOCALDEPS) $(GLOBALDEPS)
	$(FPP) $(FPPFLAGS)  -D__CHARMM -D__CHARMM_VER=$(CHARMM_MAJOR_VERSION) $(CHMSG) $(GLOBALMACROS) $(LOCALMACROS) -P $< > $@
%_tcl.c:%.i
	$(SWIGC) $(SFLAGS) -tcl -o $@ $<
%_python.c:%.i
	$(SWIGC) $(SFLAGS) -python -o $@ $<
########################## TEST PROGRAM ######################
test.F90:test.ftn $(LOCALDEPS) $(GLOBALDEPS) $(ROOT)/test.msg
	$(FPP) $(FPPFLAGS) $(DEFS) $(LOCALDEFS) -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(TESTMSG) $(DEBUG) -P $< > $@
test :: $(TESTOBJS) test.o force_look
	$(LD) $(LDFLAGS) test.o $(TESTOBJS) $(LOBJS) $(LIBS) $(MPI_LIBS) $(DEBUG) -o test
	@echo "==================================================="
	@echo "Running test :"
	@echo "==================================================="
	time ./test $(TESTARGS)	> test.out 2>&1 ; cat test.out
	@echo "==================================================="
	@echo
######################### OTHER RULES ########################
source: $(SOURCES)
$(OBJS):$(addprefix $(ROOT)/, $(addsuffix .stamp, $(DEPS)))

####################### EXECUTABLE ######################
$(EXE): $(OBJS) $(addprefix $(ROOT)/, $(addsuffix .stamp, $(DEPS)))
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

#.SILENT : %.stamp $(DISTRIB) force_look
#.SILENT : clean
clean::
	rm -f *.mod *.a *.so $(ROOT)/$(NAME).stamp $(CHARMM_SOURCES) $(SOURCES) $(OBJS) $(COBJS) $(CBINDOBJS) $(SBINDOBJS) $(DBINDOBJS) $(EXE) $(HEADS) test.F90 test.dat test.out test.o 
	if [ -f test ] ; then rm -f test ; fi
	pushd ../ >/dev/null; rm -f $(ROOT)/$(NAME).stamp ; popd >/dev/null
charmm::$(addsuffix .stamp, $(CHMDEPS))
