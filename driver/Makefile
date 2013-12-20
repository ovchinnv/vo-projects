ROOT=..
OBJS=driver.o
NAME=driver
DEPS=parser multicom output parser files timer continua dynamol
EXE=muphen

# take all object found in object directory
ALLOBJS=$(notdir $(wildcard $(ROOT)/obj/*.o) )
#exclude those that has unresolved functions or that have their own main functions
EXCOBJS=molsim.o $(OBJS)
LINKOBJS=$(filter-out $(EXCOBJS),$(ALLOBJS))

LOCALMACROS=-include codes.def
LOCALDEPS=codes.def
LOCALDEFS+=-DSIZE=gridsize

include $(ROOT)/Makefile.inc

driver : $(EXE)
	
exe : driver
	