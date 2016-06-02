LIB=libtemp_acemd.so
LIBS+=$(LIB)
LINKOBJS=\
 constants.o\
 output.o \
 parser.o \
 cmd.o prm.o \
 files.o \
 rng.o clcg.o \
 mpi_compat_stub.o \
 tempering_plugin_acemd.o \
 tempering_acemd.o \
 tempering.o 

default : acemd
acemd : $(LIB)
$(LIB) : $(LOBJS) MakeACEMD.inc

