ROOT=..
COMMON=.
OBJS=gridsize.o state.o bc.o grid.o
NAME=chest
SUBDIRS=io objects
DEPS=output parser constants files timer

LOCALDEPS=PBmacro1.def PBmacro2.def message.def

include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done

bc.o grid.o: gridsize.o state.o
state.o: gridsize.o
grid.F90: message.def

multigrid.o: grid.o state.o bc.o timer.o $(DIRS)
multigrid.F90: multigrid.src mgutil/multigrid_bc.src mgutil/mg_macro1.src mgutil/mg_resout.src \
               mgutil/GS.src mgutil/Jacobi.src mgutil/residual.src mgutil/coarsen.src mgutil/refine.src
