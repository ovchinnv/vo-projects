DIRS=bestfit constants dynamol erf files lu multicom multidiag output parselist parser random string timer vectors
STAMPS=$(addsuffix .stamp $(DIRS))

all :
	mkdir -p include
	mkdir -p obj
	mkdir -p libs
	for d in $(DIRS); do (cd $$d; $(MAKE)); done
clean ::
	for d in $(DIRS); do (cd $$d; $(MAKE) clean); done
	if [ -d include ] ; then  (cd include ; rm -f *.mod *.h ) ;  rmdir include 2>/dev/null ; fi
	if [ -d obj ] ; then  (cd obj ; rm -f *.o ) ;  rmdir obj 2>/dev/null ; fi
	if [ -d lib ] ; then  (cd lib ; rm -f *.a ) ;  rmdir lib 2>/dev/null ; fi

%:%.stamp
%.stamp:
	(cd $< ; $(MAKE))
