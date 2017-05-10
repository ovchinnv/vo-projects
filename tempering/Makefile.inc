LIB=libadtemp
LIBS+=$(LIB).so
LINKOBJS=\
 constants.o\
 output.o \
 parser.o \
 cmd.o prm.o \
 files.o \
 rng.o clcg.o \
 mpi_compat_stub.o \
 ivector.o \
 isort.o

default : plugin tcl python
plugin : $(LIB).so
$(LIB).so : $(OBJS) $(COBJS) $(LOBJS) Makefile.inc

tcl : $(LIB)_tcl.so
$(LIB)_tcl.so : tempering_tcl.o $(OBJS)
	$(LD) $(LDFLAGS) -shared tempering_tcl.o $(OBJS) $(LOBJS) $(LTCL) -Wl,-z,defs -o $@

python : _tempering.so tempering_wrap.py
_tempering.so : tempering_python.o $(OBJS)
	$(LD) $(LDFLAGS) -shared tempering_python.o $(OBJS) $(LOBJS) $(LPYTHON) -Wl,-z,defs -o $@
	cat tempering_wrap.py >> tempering.py

pythoninstall : python
	cp -f _tempering.so tempering.py $(PYTHONINSTALL)

install : pythoninstall

tempering_python.o : tempering.i
tempering_tcl.o : tempering.i
tempering.i : $(HEADS)

clean ::
	rm -f tempering_tcl.* tempering_python.* tempering.py*
