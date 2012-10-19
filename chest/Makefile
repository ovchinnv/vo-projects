ROOT=..
COMMON=.
OBJS=gridsize.o state.o bc.o grid.o multigrid.o PBmain.o ches.o
NAME=chest
SUBDIRS=io objects
DEPS=output parser constants files timer
EXE=ches

LOCALDEPS=PBmacro1.def PBmacro2.def message.def

# additional objects for linking ches executable

LINKOBJS=output.o files.o parser.o constants.o timer.o\
     plot3Dio.o chest.o util.o formats.o\
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
PBmain.o: state.o gridsize.o bc.o grid.o multigrid.o
PBmain.F90: PBmacro1.def PBmacro2.def 
multigrid.o: grid.o state.o bc.o
multigrid.F90: mgutil/multigrid_bc.src mgutil/mg_macro1.src mgutil/mg_resout.src \
               mgutil/GS.src mgutil/Jacobi.src mgutil/residual.src mgutil/coarsen.src mgutil/refine.src
grid.F90: message.def
state.o: gridsize.o
bc.o grid.o: gridsize.o state.o

clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done

test:: force_look
	$(MAKE) -C test

.SILENT : test
