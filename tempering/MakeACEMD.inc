LIB=libtemp_acemd.so
LIBS+=$(LIB)
LINKOBJS=\
 tempering_plugin_acemd.o

default : acemd
acemd : $(LIB)
$(LIB) : $(LOBJS) MakeACEMD.inc

