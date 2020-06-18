ROOT=../..
COMMON=..
OBJS=multigrid.o ches.o mgutil.o
CUOBJS=CUmgutil.o
HEADS=multigrid_c.h
NAME=continua/chest
SUBDIRS=
DEPS=output parser constants files timer multicom vectors continua/objects continua/io continua/common

LOCALDEPS=macro1.def macro2.def setcomm.def
LOCALDEPS+=source.defs
LOCALMACROS+=-include source.defs
LOCALDEFS+=-D__BINDC 
CUINCLUDES+=-I mgutil

include $(COMMON)/Makefile.inc
include $(ROOT)/Makefile.inc

ifdef __CUDA
OBJS+=$(CUOBJS)
LOCALDEFS+=-D__CUDA
endif

chest : $(EXE)

ches.o :  multigrid.o
ifdef __CUDA
multigrid.o : mgutil.o CUmgutil.o $(HEADS)
else
CUmgutil.o : mgutil/mgkernels.cu mgutil/bckernels.cu mgutil/residual.cu mgutil/coarsen.cu mgutil/refine.cu
multigrid.o : mgutil.o $(HEADS)
endif
#
# 6/20 : in GNU 10, need a -fallow-argument-mismatch to pass array sections in F77 style; the new standard way would
# require an explicit interface
multigrid.o:multigrid.F90 $(LIBDEPS)
	@echo " "$<" ==> "$@" ..."
	$(FC) $(FFLAGS) -fallow-argument-mismatch $(DEFS) $(PIC) $(LOCALDEFS) $(INCLUDES) $(MPI_INCLUDE) $< -o $@

multigrid.F90: mgutil/multigrid_bc.src $(ROOT)/parser/parser_macro1.src mgutil/mg_resout.src \
               mgutil/GS.src mgutil/Jacobi.src mgutil/residual.src mgutil/coarsen.src mgutil/refine.src\
               state/bcutil/bc_util.src
ches.F90 : $(COMMON)/outmacro1.def
clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done

test:: force_look
	$(MAKE) -C test

.SILENT : test
