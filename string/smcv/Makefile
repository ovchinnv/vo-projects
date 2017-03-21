ROOT=../..
CVOBJ=cv_posi_com.o cv_dist_com.o cv_angle_com.o cv_anglvec.o cv_dihe_com.o cv_rmsd.o cv_drmsd.o cv_proj.o cv_quaternion.o cv_qcomp.o
OBJS=cv_types.o cv_common.o cv_frames.o cv_cvrms.o smcv_master.o smcv.o smcv_add.o smcv_check_set.o smcv_stats.o $(CVOBJ)
HEADS=smcv_interface.h
NAME=string/smcv
DEPS=output vectors bestfit constants multicom random lu multidiag parser dynamol parser mpi_stub

COMMON=..

LOCALMACROS=-include source.defs -include $(COMMON)/source.defs
LOCALMACROS+=-include $(ROOT)/dmol.defs
LOCALDEPS=source.defs $(ROOT)/mpitype.def $(COMMON)/source.defs
LOCALDEPS+=$(ROOT)/dmol.defs cv_frames_calc.def
LOCALDEPS+=$(ROOT)/dmol.defs $(ROOT)/dmol_select_macro.def

include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

smcv_add.F90: $(ROOT)/dmol_select_macro.def

$(SOURCES) : $(ROOT)/source.defs

cv_frames.o : cv_common.o
cv_common.o smcv_add.o smcv.o : $(COMMON)/sm_config.mod $(COMMON)/sm_var.mod $(COMMON)/multicom_serial_string.mod
cv_common.o : cv_types.o

$(CVOBJ) cv_cvrms.o : cv_types.o cv_common.o
cv_posi_com.o cv_quaternion.o: cv_frames.o 
cv_proj.o : cv_drmsd.o
cv_qcomp.o : cv_quaternion.o
cv_cvrms.o : $(CVOBJ)
smcv.o : smcv_master.o cv_frames.o cv_quaternion.o cv_common.o smcv_stats.o smcv_add.o $(COMMON)/ftsm/ftsm_var.mod $(HEAD)
smcv_add.o : $(CVOBJ) cv_quaternion.o cv_types.o cv_common.o cv_frames.o cv_cvrms.o 
smcv.o smcv_add.o : smcv_check_set.o
smcv_master.o : cv_common.o $(CVOBJ) cv_cvrms.o smcv_stats.o

smcv.src smcv_stats.src : smcv_interface.ch
smcv.F90 smcv_stats.F90 : smcv_interface.h
smcv.F90 : $(ROOT)/exec.def $(ROOT)/contexts.def 
smcv_check_set.F90 smcv_check_set.src : smcv_check_set.def


clean::
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done
	rm -f smcv_interface.ch
