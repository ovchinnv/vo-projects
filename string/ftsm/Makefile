ROOT=../..
OBJS=ftsm_var.o ftsm_rex.o ftsm.o ftsm_voronoi.o ftsm_compute.o ftsmv2_compute.o ftsm_stats.o ftsm_io.o ftsm_util.o\
     ftsm_connect.o ftsm_min.o ftsm_rep.o ftsm_addatoms.o ftsm_interpolate.o
HEAD=interface.h
NAME=string/ftsm
COMMON=..
DEPS=constants bestfit output vectors dynamol multicom parser lu mpi_stub
INCLUDES= -I $(COMMON)

LOCALMACROS=-include source.defs -include $(COMMON)/source.defs 
LOCALMACROS+=-include $(ROOT)/dmol.defs
LOCALDEPS=source.defs $(ROOT)/mpitype.def $(COMMON)/source.defs ftsm_load_fcor.def ftsm_load_ocor.def
LOCALDEPS+= $(ROOT)/dmol.defs 
include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

ftsm_var.o: $(COMMON)/sm_config.mod $(COMMON)/sm_var.mod $(COMMON)/multicom_serial_string.mod
ftsm.o: ftsm_var.o ftsm_rex.o ftsm_voronoi.o ftsm_compute.o ftsmv2_compute.o ftsm_io.o ftsm_stats.o ftsm_util.o ftsm_min.o ftsm_rep.o\
 ftsm_connect.o ftsm_addatoms.o ftsm_interpolate.o
ftsm_connect.o: $(COMMON)/tsp.mod
ftsm_compute.o ftsmv2_compute.o ftsm_util.o ftsm_min.o ftsm_rex.o ftsm_voronoi.o ftsm_rep.o ftsm_connect.o ftsm_addatoms.o ftsm_io.o ftsm_stats.o : \
 ftsm_var.o
ftsm_stats.o: ftsm_rex.o ftsm_voronoi.o ftsm_io.o ftsm_connect.o ftsm_util.o ftsm_compute.o
ftsm_addatoms.o ftsm_io.o ftsm_min.o ftsm_stats.o ftsm_rep.o ftsm_interpolate.o : ftsm_util.o
ftsm_util.o: ftsm_compute.o ftsmv2_compute.o
ftsm_interpolate: ftsm_io.o

ftsm.src: interface.ch
ftsm.F90: interface.h

clean::
	rm -f interface.ch
