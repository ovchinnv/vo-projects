ifeq ($(WATERSHELL),1)

DEPS+=watershell
include $(ROOT)/watershell/linkobjs.inc
LINKOBJS+=watershell.o


LOCALMACROS+=-D__WATERSHELL -include watershell.defs
LOCALDEPS+=watershell.defs
endif
