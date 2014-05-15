# THIS MAKEFILE IS GENERATED AUTOMATICALLY FROM A MASTER SOURCE CODE TREE
ROOT=..
COMMON=$(ROOT)/scommon
OBJS=$(VECTORS) $(VECTORLISTS) $(SORT)
VECTORS=ivector.o rvector.o
VECTORLISTS=ivector_list.o rvector_list.o
SORT=isort.o rsort.o
LOCALDEFS=-D'__CHARMM_ONLY=!'
CHARMM_SOURCES=$(OBJS:%.o=%.src)
include $(ROOT)/Makefile.inc
#### RULES:####
charmm : $(CHARMM_SOURCES)
%.src:%.tmp $(GLOBALDEPS)
	$(FPP) $(FPPFLAGS) -D__CHARMM $(GLOBALMACROS) $< > $@
ivector.tmp: vector_template.ftn
	$(FPP) $(FPPFLAGS) -D'__DATANAME=i' -D'__DATATYPE=int' -D'__TYPENAME=int' $(CHMSG) -D__CHARMM $< >$@
ivector_list.tmp: vector_list_template.ftn
	$(FPP) $(FPPFLAGS) -P -D'__DATANAME=i' -D'__DATATYPE=int' -D'__TYPENAME=int' $(CHMSG) -D__CHARMM $< >$@
isort.tmp: sort_template.ftn
	$(FPP) $(FPPFLAGS) -P -D'__DATANAME=I' -D'__DATATYPE=int' -D'__TYPENAME=int' $(CHMSG) -D__CHARMM $< >$@

ivector.tmp: vector_template.ftn
	$(FPP) $(FPPFLAGS) -D'__DATANAME=i' -D'__DATATYPE=int' -D'__TYPENAME=int' $(CHMSG) -D__CHARMM $< >$@
ivector_list.tmp: vector_list_template.ftn
	$(FPP) $(FPPFLAGS) -P -D'__DATANAME=i' -D'__DATATYPE=int' -D'__TYPENAME=int' $(CHMSG) -D__CHARMM $< >$@
isort.tmp: sort_template.ftn
	$(FPP) $(FPPFLAGS) -P -D'__DATANAME=I' -D'__DATATYPE=int' -D'__TYPENAME=int' $(CHMSG) -D__CHARMM $< >$@
rvector.tmp: vector_template.ftn
	$(FPP) $(FPPFLAGS) -D'__DATANAME=r' -D'__DATATYPE=float' -D'__TYPENAME=real' $(CHMSG) -D__CHARMM $< >$@
rvector_list.tmp: vector_list_template.ftn
	$(FPP) $(FPPFLAGS) -P -D'__DATANAME=r' -D'__DATATYPE=float' -D'__TYPENAME=real' $(CHMSG) -D__CHARMM $< >$@
rsort.tmp: sort_template.ftn
	$(FPP) $(FPPFLAGS) -P -D'__DATANAME=R' -D'__DATATYPE=float' -D'__TYPENAME=real' $(CHMSG) -D__CHARMM $< >$@

clean ::
	rm -f $(CHARMM_SOURCES) *tmp
