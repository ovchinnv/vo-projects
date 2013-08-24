ROOT=..
VECTORS=ivector.o rvector.o
VECTORLISTS=ivector_list.o rvector_list.o
NAME=vectors
OBJS=$(VECTORS) $(VECTORLISTS)

include $(ROOT)/Makefile.inc

LOCALDEFS=-D'__CHARMM_ONLY=!'
CHARMM_SOURCES=ivector.src ivector_list.src rvector.src rvector_list.src 

#####################################################
ivector.o: ivector.F90
ivector.F90: vector_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) -U int -D'__DATANAME=i' -D'__DATATYPE=int' -D'__TYPENAME=int' $(MSG) $< >$@
rvector.o: rvector.F90
rvector.F90: vector_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) -D'__DATANAME=r' -D'__DATATYPE=float' -D'__TYPENAME=real' $(MSG) $< >$@
####################################################
ivector_list.o: ivector.o ivector_list.F90
ivector_list.F90: vector_list_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) -U int -P -D'__DATANAME=i' -D'__DATATYPE=int' -D'__TYPENAME=int' $(MSG) $< >$@
rvector_list.o: rvector.o rvector_list.F90
rvector_list.F90: vector_list_template.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) -P -D'__DATANAME=r' -D'__DATATYPE=float' -D'__TYPENAME=real' $(MSG) $< >$@
####################################################
clean::
	rm -f *.mod *.o *.a ivector.F90 rvector.F90 ivector_list.F90 rvector_list.F90 $(CHARMM_SOURCES) *.tmp \
	$(ROOT)/$(NAME).stamp

charmm: $(CHARMM_SOURCES)
%.src:%.tmp $(GLOBALDEPS)
	$(FPP) $(FPPFLAGS) -D__CHARMM $(GLOBALMACROS) $< > $@
ivector.tmp: vector_template.ftn
	$(FPP) $(FPPFLAGS) -D'__DATANAME=i' -D'__DATATYPE=int' -D'__TYPENAME=int' $(MSG) -D__CHARMM $< >$@
ivector_list.tmp: vector_list_template.ftn
	$(FPP) $(FPPFLAGS) -P -D'__DATANAME=i' -D'__DATATYPE=int' -D'__TYPENAME=int' $(MSG) -D__CHARMM $< >$@
rvector.tmp: vector_template.ftn
	$(FPP) $(FPPFLAGS) -D'__DATANAME=r' -D'__DATATYPE=float' -D'__TYPENAME=real' $(MSG) -D__CHARMM $< >$@
rvector_list.tmp: vector_list_template.ftn
	$(FPP) $(FPPFLAGS) -P -D'__DATANAME=r' -D'__DATATYPE=float' -D'__TYPENAME=real' $(MSG) -D__CHARMM $< >$@
