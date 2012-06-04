ROOT=..
OBJS=sm_var.o sm_config.o splines.o sm_util.o sm_main.o
NAME=string
SUBDIRS=smok smcv ftsm
DEPS=multicom multidiag
DEPS+=$(addprefix $(NAME)/,$(SUBDIRS))

LOCALMACROS=-include source.defs
LOCALDEPS=mpitype.def source.defs
INCLUDES=-I $(ROOT)/modules


include $(ROOT)/Makefile.inc

sm_util.o:multicom.stamp multidiag.stamp smcv

smcv:
	$(MAKE) -C ./smcv cv_common.o smcv_master.o cv_frames.o

charmm: $(CHARMM_SOURCES)

.PHONY : smcv
