ROOT=..
OBJS=sm_var.o sm_config.o splines.o sm_util.o sm_main.o tsp.o multicom_serial_string.o
NAME=string
SUBDIRS=sm0k ftsm smcv
DEPS=multicom multidiag output mpi_stub
DEPS+=$(addprefix $(NAME)/,$(SUBDIRS))
CHARMM_DIR=stringm

LOCALMACROS+=-include source.defs
LOCALDEPS+=$(ROOT)/mpitype.def source.defs 
LOCALDEPS+=dmol_select_macro.def
LOCALMACROS+=-include acemd.defs
LOCALDEPS+=acemd.defs 

include $(ROOT)/Makefile.inc

sm_util.o: splines.o multicom_serial_string.o

charmm:: $(CHARMM_SOURCES)
ifneq ($(CHARMM_DIR),.)
	mkdir -p $(CHARMM_DIR)
	/bin/sh -c 'for file in `ls *.src 2>/dev/null`; do  ln -fs `pwd`/$$file $(CHARMM_DIR)/$$file; done'
endif
	for d in $(SUBDIRS); do \
	(cd $$d; $(MAKE) charmm; \
	/bin/sh -c 'for file in `ls *.src 2>/dev/null`; do  ln -fs `pwd`/$$file ../$(CHARMM_DIR)/$$file; done'); \
	done

clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done
	rm -f $(CHARMM_DIR)/*.src 2>/dev/null
ifneq ($(CHARMM_DIR),.)
	rm -fr $(CHARMM_DIR) 2>/dev/null
endif
