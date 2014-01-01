ROOT=..
OBJS=stats.o sysinfo.o system.o sysmanip.o verlet.o dmol.o
NAME=dynamol
INCLUDES=-I $(ROOT)/include
DEPS=output parser vectors files multicom random bestfit constants datatypes
SUBDIRS=param struc coor compute
DEPS+=$(addprefix $(NAME)/,$(SUBDIRS))

# additional objects for linking a static library
LINKOBJS=constants.o output.o files.o parser.o multicom_aux.o            \
     ivector.o rng.o clcg.o stats.o system.o sysinfo.o sysmanip.o bestfit.o \
     ch_param.o atompar.o anglpar.o bondpar.o dihepar.o                     \
     psf.o psfatom.o tlist.o                                                \
     charmmio.o pdbio.o freeio.o                                            \
     angles.o bonds.o dihedrals.o


LOCALMACROS=-include source.defs
LOCALDEPS=source.defs multicom_def.src

include $(ROOT)/Makefile.inc

sysinfo.o sysmanip.o : system.o
system.o: stats.o multicom_def.src $(DIRS)
verlet.o: system.o
system.F90: parser/parser_macro1.src

libraries : libdynamol.a

clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done
