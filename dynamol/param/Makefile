ROOT=../..

OBJS = ch_param.o atompar.o bondpar.o anglpar.o dihepar.o

INCLUDES=-I $(ROOT)/modules

include $(ROOT)/Makefile.inc

ch_param.o: atompar.o bondpar.o anglpar.o dihepar.o parser output constants vectors
