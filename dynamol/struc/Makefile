ROOT=../..
OBJS=psf.o tlist.o psfatom.o
include $(ROOT)/Makefile.inc
INCLUDES=-I $(ROOT)/modules

#all: libstruc.a

#libstruc.a: psf.o
#psf.o: psf.F90 tlist.o psfatom.o  $(ROOT)/parser.o $(ROOT)/output.o $(ROOT)/size.o
tlist.o psfatom.o psf.o: parser output
psf.o : tlist.o psfatom.o

output:
	$(MAKE) -C $(ROOT)/output modules
parser:
	$(MAKE) -C $(ROOT)/parser modules

.PHONY : parser output
