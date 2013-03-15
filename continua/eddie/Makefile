ROOT=../..
COMMON=..
OBJS=eddie.o
NAME=continua/eddie
SUBDIRS=finitediff integrator pressure statistics
DEPS=output parser constants files timer multicom vectors continua/objects continua/io continua/common continua/chest
EXE=eddie

LOCALDEPS=

# additional objects for linking executable

LINKOBJS=output.o files.o parser.o constants.o timer.o\
         multicom.o multicom_aux.o \
         gridsize.o state.o grid.o datatypes.o fieldbcdef.o\
         plot3Dio.o chest.o util.o formats.o\

include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

chest : $(EXE)

ches.o: PBmain.o

PBmain.o :  multigrid.o
PBmain.F90: PBmacro1.def PBmacro2.def
multigrid.F90: mgutil/multigrid_bc.src mgutil/mg_macro1.src mgutil/mg_resout.src \
               mgutil/GS.src mgutil/Jacobi.src mgutil/residual.src mgutil/coarsen.src mgutil/refine.src
clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done

test:: force_look
	$(MAKE) -C test

.SILENT : test
