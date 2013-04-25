ROOT=../..
COMMON=..
OBJS=multigrid.o ches.o
NAME=continua/chest
SUBDIRS=
DEPS=output parser constants files timer multicom vectors continua/objects continua/io continua/common

LOCALDEPS=macro1.def macro2.def

include $(COMMON)/Makefile.inc
include $(ROOT)/Makefile.inc

chest : $(EXE)

ches.o :  multigrid.o
multigrid.F90: mgutil/multigrid_bc.src $(ROOT)/parser/parser_macro1.src mgutil/mg_resout.src \
               mgutil/GS.src mgutil/Jacobi.src mgutil/residual.src mgutil/coarsen.src mgutil/refine.src
clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done

test:: force_look
	$(MAKE) -C test

.SILENT : test
