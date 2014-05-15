# THIS MAKEFILE IS GENERATED AUTOMATICALLY FROM A MASTER SOURCE CODE TREE
ROOT=..
COMMON=$(ROOT)/scommon
OBJS=ftsm_var.o ftsm_rex.o ftsm.o ftsm_voronoi.o ftsm_compute.o ftsmv2_compute.o ftsm_stats.o ftsm_io.o ftsm_util.o\
     ftsm_connect.o ftsm_min.o ftsm_rep.o ftsm_addatoms.o
LOCALMACROS=-include source.defs -include $(COMMON)/source.defs 
LOCALDEPS=source.defs $(ROOT)/mpitype.def $(COMMON)/source.defs ftsm_load_fcor.def ftsm_load_ocor.def
CHARMM_SOURCES=$(OBJS:%.o=%.src)
include $(ROOT)/Makefile.inc
#### RULES:####
charmm : $(CHARMM_SOURCES)


clean ::
	rm -f $(CHARMM_SOURCES) *tmp
