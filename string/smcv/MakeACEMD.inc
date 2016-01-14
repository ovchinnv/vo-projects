SUBDIRS=acemd
DEPS+=$(addprefix $(NAME)/,$(SUBDIRS))
LOCALDEFS+=-D__ACEMD
