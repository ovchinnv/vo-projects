LIB=libadtemp.so
LIBS+=$(LIB)
LINKOBJS=\
 constants.o\
 output.o \
 parser.o \
 cmd.o prm.o \
 files.o \
 rng.o clcg.o \
 mpi_compat_stub.o \
 ivector.o \
 isort.o \
 tempering_plugin_acemd.o \
 tempering_plugin.o \
 tempering.o \
 energy_grid.o \
 griddata.o

default : plugin
plugin : $(LIB)
$(LIB) : $(LOBJS) Makefile.inc

