ROOT=../..
OBJS = ch_param.o atompar.o bondpar.o anglpar.o dihepar.o cmappar.o
NAME=dynamol/param
DEPS=parser output constants vectors
include $(ROOT)/Makefile.inc

ch_param.o: atompar.o bondpar.o anglpar.o dihepar.o cmappar.o
ch_param.F90: ptbl.def
