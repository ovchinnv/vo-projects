ROOT=..
OBJS=driver.o
NAME=driver
DEPS=parser multicom output parser files timer
EXE=muphen

# take all object found in object directory
ALLOBJS=$(notdir $(wildcard $(ROOT)/obj/*.o) )
#exclude those that has unresolved functions or that have their own main functions
EXCOBJS=ches.o molsim.o $(OBJS)
LINKOBJS=$(filter-out $(EXCOBJS),$(ALLOBJS))

#output.o files.o parser.o constants.o timer.o\
#     multicom.o multicom_aux.o parselist.o\
#     plot3Dio.o chest.o util.o formats.o\
#     gridsize.o state.o grid.o datatypes.o fieldbcdef.o\
#     object.o \
#     molecule.o PDB.o\
#     erf_ext.o erfappx.o erfsun_value.o erfsun_ref.o\
#     ivector_list.o \
#     ivector.o rng.o clcg.o \
#     stats.o system.o sysinfo.o sysmanip.o bestfit.o   \
#     ch_param.o atompar.o anglpar.o bondpar.o dihepar.o\
#     psf.o psfatom.o tlist.o                           \
#     mol_formats.o charmmio.o pdbio.o freeio.o         \
#     angles.o bonds.o dihedrals.o



include $(ROOT)/Makefile.inc

driver : $(EXE)
	
exe : driver
	