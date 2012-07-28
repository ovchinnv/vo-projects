
DEFS+=-D SIZE=gridsize \
      -D __TRUNCATESUMS  \
      -D __erf=erfo7 \
      -D __TIMER
#      -D __VECCOPY
#      -D __PVEC \
#      -D __erf=erfsun \
#      -D __CVEC \
#      -D __RVEC \
#      -D __erf=erfo5\
#      -D __OPENMP
# ! note: the VEC routines are probably slower
#     -D'__erf(__A)=erfsun(DBLE(__A))' # definition of Sun error function for single precision compilation

LOCALMACROS+=-include $(COMMON)/source.defs
LOCALDEPS+=$(COMMON)/source.defs
INCLUDES+=-I $(COMMON)

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
