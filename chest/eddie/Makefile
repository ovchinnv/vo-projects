ROOT=../..
COMMON=..
OBJS=plot3Dio.o util.o formats.o chest.o
NAME=chest/io
DEPS=parser output files

LOCALMACROS=-include source.defs
LOCALDEPS=source.defs

include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

plot3Dio.o chest.o : $(COMMON)/gridsize.mod

libs : libio.a
