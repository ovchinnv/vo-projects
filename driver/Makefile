ROOT=..
OBJS=driver.o
NAME=driver
DEPS=parser multicom output files timer continua continua/common continua/chest continua/io continua/objects dynamol state grid
EXE=muphen

# take all object found in object directory
ALLOBJS=$(notdir $(wildcard $(ROOT)/obj/*.o) )
#exclude those that have unresolved functions or that have their own main functions, or functions with identical/conflicting names
EXCOBJS=molsim.o watershell_plugin_acemd.o sm_plugin_acemd.o elec.o $(OBJS)
LINKOBJS=$(filter-out $(EXCOBJS),$(ALLOBJS))

LOCALMACROS=-include codes.def
LOCALDEPS=codes.def
LOCALDEFS+=-DSIZE=gridsize

include $(ROOT)/Makefile.inc

driver : $(EXE)
	
exe : driver
	