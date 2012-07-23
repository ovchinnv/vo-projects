ROOT=..
OBJS=sm_var.o sm_config.o splines.o sm_util.o sm_main.o
NAME=string
SUBDIRS=sm0k smcv ftsm
DEPS=multicom multidiag
DEPS+=$(addprefix $(NAME)/,$(SUBDIRS))
CHARMM_DIR=stringm

LOCALMACROS=-include source.defs
LOCALDEPS=mpitype.def source.defs dmol_select_macro.def

include $(ROOT)/Makefile.inc

sm_util.o : splines.o

charmm: $(CHARMM_SOURCES)
	mkdir -p $(CHARMM_DIR)
	/bin/sh -c 'for file in `ls *.src 2>/dev/null`; do  ln -fs `pwd`/$$file $(CHARMM_DIR)/$$file; done'
	for d in $(SUBDIRS); do \
	(cd $$d; $(MAKE) charmm; \
	/bin/sh -c 'for file in `ls *.src 2>/dev/null`; do  ln -fs `pwd`/$$file ../$(CHARMM_DIR)/$$file; done'); \
	done

clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done
	rm -fr $(CHARMM_DIR)
