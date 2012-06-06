ROOT=..
OBJS=sm_var.o sm_config.o splines.o sm_util.o sm_main.o
NAME=string
SUBDIRS=sm0k smcv ftsm
DEPS=multicom multidiag
DEPS+=$(addprefix $(NAME)/,$(SUBDIRS))

LOCALMACROS=-include source.defs
LOCALDEPS=mpitype.def source.defs dmol_select_macro.def

include $(ROOT)/Makefile.inc

sm_util.o : splines.o

charmm: $(CHARMM_SOURCES)

clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done

