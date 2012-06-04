ROOT=../..
OBJS=angles.o bonds.o dihedrals.o
NAME=dynamol/compute
DEPS=constants dynamol/struc dynamol/param
INCLUDES=-I $(ROOT)/modules
include $(ROOT)/Makefile.inc

$(ROOT)/$(NAME).stamp :: interface.h
	file=$< ; ln -fs `pwd`/$$file ${ROOT}/include/$$file;
	echo -n > $@

.SILENT : $(ROOT)/$(NAME).stamp
