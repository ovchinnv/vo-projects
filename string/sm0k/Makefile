ROOT=../..
OBJS=sm0k.o
NAME=string/sm0k
DEPS=constants bestfit output dynamol/struc dynamol/coor multicom parser
COMMON=..
INCLUDES= -I $(COMMON) -I $(ROOT)/modules

LOCALMACROS=-include ../source.defs

include $(ROOT)/Makefile.inc

sm0k.o: $(COMMON)/sm_config.o
	if [ -f sm_config.mod ] ; then mv -f sm_config.mod ../sm_config.mod 2>/dev/null; fi

charmm: $(CHARMM_SOURCES)
