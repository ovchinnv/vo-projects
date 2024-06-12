ifneq ($(WATERSHELL),0)
include ../watershell/version.inc

ifeq ($(version),2)
 LINKOBJS+=watershell2.o
else
 version=
 LINKOBJS+=watershell.o
endif

DEPS+=watershell
include $(ROOT)/watershell/linkobjs.inc
LINKOBJS+=ugrid.o

LOCALMACROS+=-D__WATERSHELL -D__WATERSHELL_VERSION=$(version)

plugin_master.F90 : $(ROOT)/watershell.stamp  watershell_vars.src watershell_init.src watershell_done.src

watershell_vars.src : force_look
	sed -n '/!#vars 0/{:a;n;/!#end vars 0/b;p;ba}' $(ROOT)/watershell/watershell_plugin.F90 >$@

watershell_init.src : force_look
	sed -n '/!#init 0/{:a;n;/!#end init 0/b;p;ba}' $(ROOT)/watershell/watershell_plugin.F90 >$@

watershell_done.src : force_look
	sed -n '/!#done 0/{:a;n;/!#end done 0/b;p;ba}' $(ROOT)/watershell/watershell_plugin.F90 >$@

.SILENT : *src

endif
LOCALMACROS+=-include watershell.defs
LOCALDEPS+=watershell.defs
