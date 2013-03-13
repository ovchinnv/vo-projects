ROOT=../..
COMMON=..
OBJS=multigrid.o PBmain.o ches.o
NAME=continua/chest
SUBDIRS=
DEPS=output parser constants files timer multicom vectors continua/objects continua/io continua/common
EXE=ches

LOCALDEPS=PBmacro1.def PBmacro2.def
#LOCALDEPS+=$(COMMON)/grid.mod $(COMMON)/gridsize.mod $(COMMON)/datatypes.mod $(COMMON)/fieldbcdef.mod


# additional objects for linking ches executable

LINKOBJS=output.o files.o parser.o constants.o timer.o\
     multicom.o multicom_aux.o parselist.o\
     plot3Dio.o chest.o util.o formats.o\
     gridsize.o state.o grid.o datatypes.o fieldbcdef.o\
     object.o \
     molecule.o PDB.o\
     erf_ext.o erfappx.o erfsun_value.o erfsun_ref.o\
     ivector_list.o \
     ivector.o rng.o clcg.o \
     stats.o system.o sysinfo.o sysmanip.o bestfit.o   \
     ch_param.o atompar.o anglpar.o bondpar.o dihepar.o\
     psf.o psfatom.o tlist.o                           \
     mol_formats.o charmmio.o pdbio.o freeio.o         \
     angles.o bonds.o dihedrals.o

include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

default : $(EXE)

ches.o: PBmain.o

#multigrid.o PBmain.o: $(COMMON)/gridsize.mod $(COMMON)/fieldbcdef.mod $(COMMON)/datatypes.mod $(COMMON)/state.mod $(COMMON)/grid.mod
PBmain.o :  multigrid.o
PBmain.F90: PBmacro1.def PBmacro2.def
multigrid.F90: mgutil/multigrid_bc.src mgutil/mg_macro1.src mgutil/mg_resout.src \
               mgutil/GS.src mgutil/Jacobi.src mgutil/residual.src mgutil/coarsen.src mgutil/refine.src
clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done

test:: force_look
	$(MAKE) -C test

.SILENT : test
