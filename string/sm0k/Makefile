ROOT=../..
OBJS=sm0k.o
NAME=string/sm0k
DEPS=constants bestfit output dynamol multicom parser
COMMON=..

LOCALMACROS=-include ../source.defs -include $(ROOT)/dmol.defs
LOCALDEPS=source.defs ../source.defs $(ROOT)/dmol.defs

include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

sm0k.o: $(COMMON)/sm_config.mod
