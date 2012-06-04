ROOT=../..
OBJS=ftsm_var.o ftsm_rex.o ftsm.o
NAME=string/ftsm
COMMON=..
DEPS=constants vectors
INCLUDES= -I $(COMMON)

LOCALMACROS=-include source.defs -include ../source.defs
LOCALDEPS=source.defs ../source.defs

include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

ftsm_var.o: $(COMMON)/sm_config.mod $(COMMON)/sm_var.mod
ftsm.o: ftsm_var.o ftsm_rex.o
ftsm_rex.o: ftsm_var.o

