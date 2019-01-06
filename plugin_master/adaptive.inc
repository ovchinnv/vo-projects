ifeq ($(ADAPTIVE),1)
ifeq ($(STRUNA),1)

DEPS+=adaptive
LINKOBJS+=adaptive_kforce.o gridscalar_log.o gridscalartype.o

LOCALMACROS+=-D__ADAPTIVE_KCV 
LOCALMACROS+=-include adaptive.defs
LOCALDEPS+=adaptive.defs

endif
endif

