ROOT=..
OBJS=rng.o clcg.o
NAME=random
DEPS=constants output parser

include $(ROOT)/Makefile.inc

rng.o : clcg.o 
clcg.o : clcg.f clcg.h
