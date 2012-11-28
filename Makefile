DIRS=bestfit constants dynamol erf files lu multicom multidiag output parselist parser random string timer vectors chest confcons enm
DISTRIB=include obj lib

all : $(DISTRIB)
	for d in $(DIRS); do (cd $$d; $(MAKE)); done
clean ::
	for d in $(DIRS); do (cd $$d; $(MAKE) clean); done
	for d in $(DISTRIB); do if [ -d $$d ] ; then rm -rf $$d 2>/dev/null ; fi ; done

$(DIRS) :: $(DISTRIB)
	(cd $@ ; $(MAKE))
$(DISTRIB) :
	if [ ! -d $@ ] ; then mkdir $@ ; fi
charmm :
		for d in $(DIRS); do (cd $$d; $(MAKE) charmm); done
.SILENT : charmm