INCLUDES+= -I $(COMMON)

$(COMMON)/%.mod : $(COMMON)/%.o
	file=$(notdir $@) ; if [ -f $$file ] ; then mv -f $$file $(COMMON)/$$file 2>/dev/null; fi
	touch -m $@ 2>/dev/null ; # change modification to update dependency status
charmm: $(CHARMM_SOURCES)

