ROOT=..
OBJS=rng.o clcg.o
NAME=random
DEPS=constants output parser

LINKOBJS=output.o parser.o cmd.o prm.o mpi_compat_stub.o

include $(ROOT)/Makefile.inc

rng.o : clcg.o 
clcg.o : clcg.f clcg.h
