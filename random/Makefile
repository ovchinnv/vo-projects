ROOT=..
OBJS=rng.o clcg.o
NAME=random
DEPS=constants output parser

INCLUDES=-I $(ROOT)/modules

include $(ROOT)/Makefile.inc

rng.o : clcg.o 
clcg.o : clcg.f clcg.h
