ROOT=../..
OBJS=ftsm_var.o ftsm_rex.o ftsm.o ftsm_voronoi.o
NAME=string/ftsm
COMMON=..
DEPS=constants bestfit output vectors dynamol multicom parser
INCLUDES= -I $(COMMON)

LOCALMACROS=-include source.defs -include ../source.defs -include $(ROOT)/dmol.defs
LOCALDEPS=source.defs ../source.defs $(ROOT)/dmol.defs

include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

ftsm_var.o: $(COMMON)/sm_config.mod $(COMMON)/sm_var.mod
ftsm.o: ftsm_var.o ftsm_rex.o ftsm_voronoi.o
ftsm_rex.o: ftsm_var.o
ftsm_voronoi.o: ftsm_var.o
