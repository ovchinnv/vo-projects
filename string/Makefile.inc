INCLUDES+= -I $(COMMON)

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
charmm: $(CHARMM_SOURCES)

