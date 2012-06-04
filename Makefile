DIRS=bestfit constants dynamol erf files lu multicom multidiag output parselist parser random string timer vectors
DISTRIB=include obj libs
STAMPS=$(addsuffix .stamp, $(DIRS))


all : $(DISTRIB)
	for d in $(DIRS); do (cd $$d; $(MAKE)); done
clean ::
	for d in $(DIRS); do (cd $$d; $(MAKE) clean); done
	for d in $(DIST); do if [ -d $$d ] ; then rm -rf $$d 2>/dev/null ; fi ; done

$(DIRS) :: $(DISTRIB)
	(cd $@ ; $(MAKE))

$(DISTRIB) :
	if [ ! -d $@ ] ; then mkdir $@ ; fi
