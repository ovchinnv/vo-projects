ROOT=..
NAME=elec
OBJS=$(addsuffix .o, vars setup filt3 direct3 direct3_nosym util)
#OBJS=vars.o setup.o filt3.o direct3.o direct3_nosym.o util.o
DEPS=parser constants output

TESTOBJS:=$(filter-out run.o, $(OBJS)) test.o

LOCALDEPS=source.def filters.def
LOCALMACROS+=$(addprefix -include , $(LOCALDEPS))

LINKOBJS= $(NAME).o

EXE=run

include $(ROOT)/Makefile.inc

.PRECIOUS : test $(EXE)

$(OBJS) : vars.o
#test : $(TESTOBJS)
#	$(LD) $(TESTOBJS) -o test
#	./test
