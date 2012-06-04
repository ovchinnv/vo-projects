ROOT=../..
OBJS = ch_param.o atompar.o bondpar.o anglpar.o dihepar.o
INCLUDES=-I $(ROOT)/modules
NAME=dynamol/param
DEPS=parser output constants vectors
include $(ROOT)/Makefile.inc

ch_param.o: atompar.o bondpar.o anglpar.o dihepar.o
