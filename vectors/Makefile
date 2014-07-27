ROOT=..
OBJS=$(VECTORS) $(VECTORLISTS) $(SORT)
VECTORS=ivector.o rvector.o
VECTORLISTS=ivector_list.o rvector_list.o
SORT=isort.o rsort.o
NAME=vectors
DEPS=output

include $(ROOT)/Makefile.inc


LOCALDEFS=-D'__CHARMM_ONLY=!'
CHARMM_SOURCES=ivector.src ivector_list.src rvector.src rvector_list.src isort.src rsort.src
$(ROOT)/$(NAME).stamp:: $(CHARMM_SOURCES)

LINKOBJS=output.o
#####################################################
ivector.o:  isort.o ivector.F90
ivector.F90: vector_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) -U int -D'__DATANAME=i' -D'__DATATYPE=int' -D'__TYPENAME=int' $(MSG) $< >$@
rvector.o: rsort.o rvector.F90
rvector.F90: vector_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) -D'__DATANAME=r' -D'__DATATYPE=float' -D'__TYPENAME=real' $(MSG) $< >$@
####################################################
ivector_list.o: ivector.o ivector_list.F90
ivector_list.F90: vector_list_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) -U int -P -D'__DATANAME=i' -D'__DATATYPE=int' -D'__TYPENAME=int' $(MSG) $< >$@
rvector_list.o: rvector.o rvector_list.F90
rvector_list.F90: vector_list_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) -P -D'__DATANAME=r' -D'__DATATYPE=float' -D'__TYPENAME=real' $(MSG) $< >$@
#####################################################
isort.o: isort.F90
isort.F90: sort_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) $(GLOBALMACROS) -U int -D'__DATANAME=I' -D'__DATATYPE=int' -D'__TYPENAME=int' $(MSG) $< >$@
rsort.o: rsort.F90
rsort.F90: sort_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) $(GLOBALMACROS) -D'__DATANAME=R' -D'__DATATYPE=float' -D'__TYPENAME=real' $(MSG) $< >$@
####################################################
clean::
	rm -f *.mod *.o *.a ivector.F90 rvector.F90 ivector_list.F90 isort.F90 rvector_list.F90 isort.F90 rsort.F90 $(CHARMM_SOURCES) *.tmp \
	$(ROOT)/$(NAME).stamp

charmm:: $(CHARMM_SOURCES)
%.src:%.tmp $(GLOBALDEPS)
	$(FPP) $(FPPFLAGS) -D__CHARMM $(GLOBALMACROS) $< > $@
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
