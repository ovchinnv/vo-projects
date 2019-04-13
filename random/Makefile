ROOT=..
OBJS=rng.o clcg.o
NAME=random
DEPS=constants output parser

LINKOBJS=output.o parser.o cmd.o prm.o mpi_compat_stub.o
TESTOBJS=$(OBJS)

include $(ROOT)/Makefile.inc

rng.o : clcg.o 
clcg.o : clcg.F clcg.h
