ROOT=..
OBJS=pnm.o
NAME=pnm
DEPS=vectors parser output constants dynamol
CHMDEPS=vectors

LOCALMACROS=-include source.defs
LOCALMACROS+=-include $(ROOT)/dmol.defs

LOCALDEPS=source.defs
LOCALDEPS+= $(ROOT)/dmol.defs 

include $(ROOT)/Makefile.inc

charmm:: vectors $(CHARMM_SOURCES)
vectors : $(ROOT)/getcode force_look
	/bin/sh -c 'source $(ROOT)/getcode;\
	get_typedef_from_file int_vector $(ROOT)/vectors/ivector.src > vectors_var ;\
	get_typedef_from_file real_vector $(ROOT)/vectors/rvector.src >> vectors_var ;\
	get_param_from_file expand_incr $(ROOT)/vectors/rvector.src >> vectors_var ;\
	get_subroutine_from_file int_vector_init $(ROOT)/vectors/ivector.src > vectors_sub ;\
	get_subroutine_from_file int_vector_done $(ROOT)/vectors/ivector.src >> vectors_sub ;\
	get_subroutine_from_file int_vector_expand $(ROOT)/vectors/ivector.src >> vectors_sub ;\
	get_function_from_file int_vector_add $(ROOT)/vectors/ivector.src >> vectors_sub ;\
	get_subroutine_from_file real_vector_init $(ROOT)/vectors/rvector.src >> vectors_sub ;\
	get_subroutine_from_file real_vector_done $(ROOT)/vectors/rvector.src >> vectors_sub ;\
	get_subroutine_from_file real_vector_expand $(ROOT)/vectors/rvector.src >> vectors_sub ;\
	get_function_from_file real_vector_add $(ROOT)/vectors/rvector.src >> vectors_sub ;\
	'

clean ::
	rm -f *.mod *.o *.a $(CHARMM_SOURCES) vectors_var vectors_sub \
	$(ROOT)/$(NAME).stamp

.SILENT : vectors
