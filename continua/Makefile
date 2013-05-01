ROOT=..
COMMON=.
OBJS=
NAME=continua
SUBDIRS=io common objects chest # eddie
DEPS=
LOCALDEPS=

include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

$(ROOT)/$(NAME).stamp:: $(SUBDIRS)

$(SUBDIRS) : $(addprefix $(NAME)/, $(addsuffix .stamp, $(SUBDIRS)))

chest : force_look
	cd chest; $(MAKE) chest

clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done

test:: force_look
	$(MAKE) -C test



.SILENT : test
