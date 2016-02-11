DIRS=driver state bestfit constants dynamol erf files lu multicom multidiag output parselist parser random string timer vectors \
confcons enm continua pnm mpi_stub
DISTRIB=include obj lib

all : $(DIRS)
clean ::
	for d in $(DIRS); do (cd $$d; $(MAKE) clean); done
	for d in $(DISTRIB); do if [ -d $$d ] ; then rm -rf $$d 2>/dev/null ; fi ; done
	rm -f *.stamp
$(DIRS) :: $(DISTRIB)
	(cd $@ ; $(MAKE))
$(DISTRIB) :
	if [ ! -d $@ ] ; then mkdir $@ ; fi
exe :: $(DISTRIB)
	(cd driver ; $(MAKE) exe)
charmm :
		for d in $(DIRS); do (cd $$d; $(MAKE) charmm); done
.SILENT : charmm
.PHONY : $(DIRS)
include Makefile.charmm # rules to build directory for patching CHARMM v. 39
