ROOT=..
OBJS=stats.o sysinfo.o system.o sysmanip.o verlet.o
INCLUDES=-I $(ROOT)/include
DIRS=param struc coor compute

# for making a static library
LINKOBJS=constants.o output.o files.o parser.o multicom_aux.o            \
     ivector.o rng.o clcg.o stats.o system.o sysinfo.o sysmanip.o bestfit.o \
     ch_param.o atompar.o anglpar.o bondpar.o dihepar.o                     \
     psf.o psfatom.o tlist.o                                                \
     charmmio.o pdbio.o freeio.o                                            \
     angles.o bonds.o dihedrals.o

include $(ROOT)/Makefile.inc

sysinfo.o sysmanip.o : system.o output
system.o: stats.o parser output vectors files multicom_aux random bestfit multicom_def.src $(DIRS)
verlet.o: system.o parser output constants random

libraries : libdynamol.a
libdynamol.a: $(OBJS)  force_look                      #force relinking
	$(AR) $(AFLAGS) $@ $(OBJS) $(LOBJS)
	/bin/sh -c 'ln -fs `pwd`/libdynamol.a ${ROOT}/lib/libdynamol.a'
clean: 
	for d in $(DIRS); do (cd $$d; $(MAKE) clean); done
	rm -f *.mod *.a $(OBJS) $(SOURCES)
