ROOT=../..
COMMON=..
OBJS=plot3dio.o util.o formats.o chestio.o fieldio.o opendxio.o
NAME=continua/io
DEPS=parser output files mpi_stub

include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

#libs : libio.a
fieldio.o : plot3dio.o util.o formats.o chestio.o opendxio.o

