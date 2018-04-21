ROOT=..
NAME=elec
OBJS=$(addsuffix .o, vars setup filt3 direct3 direct3_nosym util $(NAME))
#OBJS=vars.o setup.o filt3.o direct3.o direct3_nosym.o util.o
DEPS=parser constants output files vectors

TESTOBJS=$(filter-out $(NAME).o, $(OBJS))
TESTARGS=input

LOCALDEPS=source.def filters.def
LOCALMACROS+=$(addprefix -include , $(LOCALDEPS))

LINKOBJS= constants.o parser.o output.o cmd.o prm.o files.o ivector.o isort.o

EXE=run


include $(ROOT)/Makefile.inc

.PRECIOUS : test $(EXE)

$(filter-out vars.o, $OBJS) : vars.o
