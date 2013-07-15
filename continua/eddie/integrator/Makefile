ROOT=../../..
COMMON=../..
EDDIE=..
#
OBJS=ABintegrator.o timestep.o baseintegrator.o integrator.o
NAME=continua/eddie/integrator
DEPS=output parser constants continua/common continua/eddie/discretization
#
include $(EDDIE)/Makefile.inc
include $(COMMON)/Makefile.inc
include $(ROOT)/Makefile.inc

ABintegrator.o : baseintegrator.o timestep.o
integrator.o : ABintegrator.o 
integrator.F90 : $(EDDIE)/init_class_macro.src
