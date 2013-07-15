ROOT=../../..
COMMON=../..
EDDIE=..
OBJS=fd2c.o
OBJS=discretization.o basediscretization.o
NAME=continua/eddie/discretization
DEPS=constants continua/common

include $(EDDIE)/Makefile.inc
include $(COMMON)/Makefile.inc
include $(ROOT)/Makefile.inc

discretization.F90 : $(EDDIE)/init_class_macro.src
discretization.o : basediscretization.o
basediscretization.F90 : parser_macro1.src