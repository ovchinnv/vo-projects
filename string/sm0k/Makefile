ROOT=../..
OBJS=sm0k.o
NAME=string/sm0k
DEPS=constants bestfit output dynamol multicom parser confcons mpi_stub
COMMON=..

LOCALMACROS=-include $(COMMON)/source.defs
LOCALMACROS+= -include $(ROOT)/dmol.defs
LOCALDEPS=$(ROOT)/mpitype.def $(COMMON)/source.defs 
LOCALDEPS+=$(ROOT)/dmol.defs

include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

sm0k.o: $(COMMON)/sm_config.mod $(ROOT)/dmol_select_macro.def $(COMMON)/multicom_serial_string.mod
