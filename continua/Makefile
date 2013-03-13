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

clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done

test:: force_look
	$(MAKE) -C test



.SILENT : test
