ROOT=../..
OBJS=sm0k.o
NAME=string/sm0k
DEPS=constants bestfit output dynamol multicom parser confcons
COMMON=..

LOCALMACROS=-include $(COMMON)/source.defs
LOCALMACROS+= -include $(ROOT)/dmol.defs
LOCALDEPS=$(ROOT)/mpitype.def $(COMMON)/source.defs 
LOCALDEFS+=$(ROOT)/dmol.defs

include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

sm0k.o: $(COMMON)/sm_config.mod
