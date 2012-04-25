ROOT=..
OBJS=sm_var.o sm_config.o splines.o sm_util.o sm_main.o

LOCALMACROS=-include source.defs
LOCALDEPS=mpitype.def source.defs

INCLUDES=-I $(ROOT)/modules

include $(ROOT)/Makefile.inc

sm_util.o:multicom_aux
multicom_aux:
	$(MAKE) -C $(ROOT)/multicom multicom_aux.o modules

charmm: $(CHARMM_SOURCES)

.PHONY : multicom_aux
