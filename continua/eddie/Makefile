ROOT=../..
COMMON=..
OBJS=eddie_grid.o eddie_info.o
NAME=continua/eddie
DEPS=output parser constants continua/common


include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc
#include $(COMMON)/eddie/Makefile.inc


eddie_grid.F90: $(ROOT)/parser/parser_macro1.src
eddie_info.o: eddie_grid.o
eddie.o: eddie_grid.o eddie_info.o
eddie.F90: $(ROOT)/parser/parser_macro1.src


test:: force_look
	$(MAKE) -C test


.SILENT : test
