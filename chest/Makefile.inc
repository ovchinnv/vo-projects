
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

$(COMMON)/%.mod : $(COMMON)/%.o
	file=$(notdir $@) ; if [ -f $$file ] ; then mv -f $$file $(COMMON)/$$file 2>/dev/null; fi
	touch -m $@ 2>/dev/null ; # change modification to update dependency status
