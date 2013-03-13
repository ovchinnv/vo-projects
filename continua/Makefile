ROOT=..
COMMON=.
NAME=continua
SUBDIRS=common io objects chest # eddie
DEPS=
LOCALDEPS=

include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

default : all
all : $(SUBDIRS)
	for d in $(SUBDIRS); do (cd $$d; $(MAKE)); done

$(SUBDIRS) ::
	(cd $@ ; $(MAKE))

clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done

test:: force_look
	$(MAKE) -C test



.SILENT : test
