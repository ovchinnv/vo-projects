ROOT=../..
OBJS=sm0k.o

COMMON=..
INCLUDES= -I $(COMMON) -I $(ROOT)/modules

LOCALMACROS=-include ../source.defs

include $(ROOT)/Makefile.inc

sm0k.o: constants bestfit output psf multicom_aux sm_common_vars parser coor

sm_common_vars:
	$(MAKE) -C $(COMMON) sm_config.o

charmm: $(CHARMM_SOURCES)


.PHONY : sm_common_vars
