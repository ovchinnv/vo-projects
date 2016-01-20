INCLUDES+= -I $(COMMON) -I $(COMMON)/acemd -I $(COMMON)/ftsm -I $(COMMON)/smcv

LOCALMACROS+=-include $(COMMON)/acemd.defs
LOCALDEPS+=$(COMMON)/acemd.defs

$(COMMON)/%.mod : $(COMMON)/%.o
	file=$(notdir $@) ; \
	dir=$(dir $@) ; \
	if [ -f $$file ] ; then mv -f $$file $@ 2>/dev/null; touch -m $@ 2>/dev/null; else \
	FILE=`echo $${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ; \
	if [ -f $$FILE ] ; then mv -f $$FILE $$dir/$$FILE 2>/dev/null; \
	touch -m $$dir/$$FILE 2>/dev/null;\
	fi ; fi
# touch used to change modification to update dependency status; 
# note: pathscale converts module names to uppercase
charmm:: $(CHARMM_SOURCES)

