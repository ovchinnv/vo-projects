ROOT=../..
OBJS=psf.o tlist.o psfatom.o
INCLUDES=-I $(ROOT)/modules
include $(ROOT)/Makefile.inc

tlist.o psfatom.o psf.o: parser output
psf.o : tlist.o psfatom.o
