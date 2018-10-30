ROOT=..
OBJS=$(VECTORS) $(VECTORLISTS) $(SORT)
VECTORS=ivector.o rvector.o i8vector.o
VECTORLISTS=ivector_list.o rvector_list.o 
#i8vector_list.o
SORT=isort.o rsort.o isort_ind.o rsort_ind.o i8sort.o
NAME=vectors
DEPS=output

include $(ROOT)/Makefile.inc


LOCALDEFS=-D'__CHARMM_ONLY=!'
CHARMM_SOURCES=ivector.src ivector_list.src rvector.src rvector_list.src i8vector.src isort.src rsort.src i8sort.src rsort_ind.src isort_ind.src 
# i8vector_list.src
#$(ROOT)/$(NAME).stamp:: $(CHARMM_SOURCES)

LINKOBJS=output.o
#####################################################
ivector.o:  isort.o ivector.F90
ivector.F90: vector_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) -U int -D'__DATANAME=i' -D'__DATATYPE=int' -D'__TYPENAME=int' $(MSG) $< >$@
rvector.o: rsort.o rvector.F90
rvector.F90: vector_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) -D'__DATANAME=r' -D'__DATATYPE=float' -D'__TYPENAME=real' $(MSG) $< >$@
i8vector.o:  i8sort.o i8vector.F90
i8vector.F90: vector_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) -U int -D'__DATANAME=i8' -D'__DATATYPE=INTEGER*8' -D'__TYPENAME=int8' $(MSG) $< >$@
####################################################
ivector_list.o: ivector.o ivector_list.F90
ivector_list.F90: vector_list_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) -U int -P -D'__DATANAME=i' -D'__DATATYPE=int' -D'__TYPENAME=int' $(MSG) $< >$@
rvector_list.o: rvector.o rvector_list.F90
rvector_list.F90: vector_list_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) -P -D'__DATANAME=r' -D'__DATATYPE=float' -D'__TYPENAME=real' $(MSG) $< >$@
#i8vector_list.o: i8vector.o i8vector_list.F90
#i8vector_list.F90: vector_list_template.ftn
#	$(FPP) $(FPPFLAGS) $(DEFS) -U int -P -D'__DATANAME=i8' -D'__DATATYPE=INTEGER*8' -D'__TYPENAME=int8' $(MSG) $< >$@
#####################################################
isort.o: isort.F90
isort.F90: sort_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) $(GLOBALMACROS) -U int -D'__DATANAME=I' -D'__DATATYPE=int' -D'__TYPENAME=int' $(MSG) $< >$@
isort_ind.o: isort_ind.F90
isort_ind.F90: sort_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) $(GLOBALMACROS) -U int -D'__DATANAME=I' -D'__DATATYPE=int' -D'__TYPENAME=int' -D__INDEX $(MSG) $< >$@
rsort.o: rsort.F90
rsort.F90: sort_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) $(GLOBALMACROS) -D'__DATANAME=R' -D'__DATATYPE=float' -D'__TYPENAME=real' $(MSG) $< >$@
rsort_ind.o: rsort_ind.F90
rsort_ind.F90: sort_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) $(GLOBALMACROS) -D'__DATANAME=R' -D'__DATATYPE=float' -D'__TYPENAME=real' -D__INDEX $(MSG) $< >$@
i8sort.o: i8sort.F90
i8sort.F90: sort_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) $(GLOBALMACROS) -U int -D'__DATANAME=I8' -D'__DATATYPE=INTEGER*8' -D'__TYPENAME=int8' $(MSG) $< >$@
####################################################
clean::
	rm -f *.mod *.o *.a ivector.F90 rvector.F90 ivector_list.F90 isort.F90 rvector_list.F90 isort.F90 rsort.F90 isort_ind.F90 rsort_ind.F90 $(CHARMM_SOURCES) *.tmp \
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
isort_ind.tmp: sort_template.ftn
	$(FPP) $(FPPFLAGS) -P -D'__DATANAME=I' -D'__DATATYPE=int' -D'__TYPENAME=int' -D__INDEX $(CHMSG) -D__CHARMM $< >$@
rvector.tmp: vector_template.ftn
	$(FPP) $(FPPFLAGS) -D'__DATANAME=r' -D'__DATATYPE=float' -D'__TYPENAME=real' $(CHMSG) -D__CHARMM $< >$@
rvector_list.tmp: vector_list_template.ftn
	$(FPP) $(FPPFLAGS) -P -D'__DATANAME=r' -D'__DATATYPE=float' -D'__TYPENAME=real' $(CHMSG) -D__CHARMM $< >$@
rsort.tmp: sort_template.ftn
	$(FPP) $(FPPFLAGS) -P -D'__DATANAME=R' -D'__DATATYPE=float' -D'__TYPENAME=real' $(CHMSG) -D__CHARMM $< >$@
rsort_ind.tmp: sort_template.ftn
	$(FPP) $(FPPFLAGS) -P -D'__DATANAME=R' -D'__DATATYPE=float' -D'__TYPENAME=real' -D__INDEX $(CHMSG) -D__CHARMM $< >$@
i8vector.tmp: vector_template.ftn
	$(FPP) $(FPPFLAGS) -D'__DATANAME=i8' -D'__DATATYPE=INTEGER*8' -D'__TYPENAME=int8' $(CHMSG) -D__CHARMM $< >$@
#i8vector_list.tmp: vector_list_template.ftn
#	$(FPP) $(FPPFLAGS) -P -D'__DATANAME=i8' -D'__DATATYPE=INTEGER*8' -D'__TYPENAME=int8' $(CHMSG) -D__CHARMM $< >$@
i8sort.tmp: sort_template.ftn
	$(FPP) $(FPPFLAGS) -P -D'__DATANAME=I8' -D'__DATATYPE=INTEGER*8' -D'__TYPENAME=int8' $(CHMSG) -D__CHARMM $< >$@
