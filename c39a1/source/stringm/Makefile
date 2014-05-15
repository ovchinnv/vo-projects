SUBDIRS=bestfit lu multicom multidiag confcons parselist scommon ftsm sm0k smcv vectors
CHARMM_DIR=.
PREFLX2CPP=./prefx2cpp
charmm:: $(CHARMM_SOURCES)
ifneq ($(CHARMM_DIR),.)
	mkdir -p $(CHARMM_DIR)
	/bin/sh -c 'for file in `ls *.src 2>/dev/null`; do  ln -fs `pwd`/$$file $(CHARMM_DIR)/$$file; done'
endif
	for d in $(SUBDIRS); do \
	(cd $$d; $(MAKE) charmm; \
	/bin/sh -c 'for file in `ls *.src 2>/dev/null`; do  ln -fs `pwd`/$$file ../$(CHARMM_DIR)/$$file; done'); \
	done


charmm::
	/bin/sh -c 'for file in `ls *.src 2>/dev/null`;\
	do \
	flag=`grep -i "automatically protect all code" $$file`;\
	if [ -z "$$flag" ];then \
	echo "##IF STRINGM ! automatically protect all code" > _tmp;\
	cat $$file >>_tmp;\
	echo "##ENDIF ! automatically protect all code" >> _tmp;\
	$(PREFLX2CPP) _tmp > $$file;\
	rm -f _tmp __tmp;\
	fi;\
	done ;'\

install : charmm
all : install

clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done
	rm -f $(CHARMM_DIR)/*.src 2>/dev/null
ifneq ($(CHARMM_DIR),.)
	rm -fr $(CHARMM_DIR) 2>/dev/null
endif

