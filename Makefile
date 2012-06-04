DIRS=bestfit constants dynamol erf files lu multicom multidiag output parselist parser random string timer vectors
DIST=include obj libs
STAMPS=$(addsuffix .stamp, $(DIRS))


all : $(DIST)
	for d in $(DIRS); do (cd $$d; $(MAKE)); done
clean ::
	for d in $(DIRS); do (cd $$d; $(MAKE) clean); done
	for d in $(DIST); do if [ -d $$d ] ; then rm -rf $$d 2>/dev/null ; fi ; done

$(DIRS) :: $(DIST)
	(cd $@ ; $(MAKE))

$(DIST) :
	if [ ! -d $@ ] ; then mkdir $@ ; fi
