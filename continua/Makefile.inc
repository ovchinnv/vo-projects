
DEFS+=-DSIZE=gridsize \
      -U__TRUNCATESUMS  \
      -D__erf=__erfo7 \
#      -D__TIMER \
#      -D__VECCOPY
#      -D__PVEC \
#      -D__erf=erfsun \
#      -D__CVEC \
#      -D__RVEC \
#      -D__erf=erfo5\
#      -D__OPENMP
# ! note: the VEC routines are probably slower
#     -D'__erf(__A)=erfsun(DBLE(__A))' # definition of Sun error function for single precision compilation

LOCALMACROS+=-include $(COMMON)/source.defs
LOCALDEPS+=$(COMMON)/source.defs
INCLUDES+=-I $(COMMON)/common

$(COMMON)/%.mod : $(COMMON)/%.o
	file=$(notdir $@) ; \
	if [ -f $$file ] ; then mv -f $$file $(COMMON)/$$file 2>/dev/null; touch -m $@ 2>/dev/null ; else \
	FILE=`echo $${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ; \
	if [ -f $$FILE ] ; then mv -f $$FILE $(COMMON)/$$FILE 2>/dev/null; \
	touch -m $(COMMON)/$$FILE 2>/dev/null ; \
	rm -f $(COMMON)/$$file ; \
	ln -s $$FILE $(COMMON)/$$file ; fi ; fi
# touch used to change modification to update dependency status; 
# note: pathscale converts module names to uppercase
