ROOT=../..

OBJS = charmmio.o freeio.o pdbio.o

INCLUDES=-I $(ROOT)/modules

include $(ROOT)/Makefile.inc

$(OBJS): parser output psf

psf:
	$(MAKE) -C $(ROOT)/dynamol/struc psf.o modules
