ROOT=../..
CVOBJ=cv_posi_com.o cv_dist_com.o cv_angle_com.o cv_anglvec.o cv_dihe_com.o cv_rmsd.o cv_drmsd.o cv_proj.o cv_quaternion.o cv_qcomp.o cv_cvrms.o
OBJS=cv_types.o cv_common.o cv_frames.o cv_cvrms.o smcv_master.o smcv.o smcv_add.o $(CVOBJ)
NAME=string/smcv
DEPS=output vectors bestfit constants multicom random lu multidiag parser dynamol

COMMON=..

LOCALMACROS=-include source.defs -include $(COMMON)/source.defs 
LOCALMACROS+=-include $(ROOT)/dmol.defs
LOCALDEPS=source.defs $(ROOT)/mpitype.def $(COMMON)/source.defs 
LOCALDEPS+=$(ROOT)/dmol.defs

include $(ROOT)/Makefile.inc
include $(COMMON)/Makefile.inc

cv_frames.o : cv_common.o
cv_common.o smcv_add.o smcv.o : $(COMMON)/sm_config.mod $(COMMON)/sm_var.mod $(COMMON)/mpi_compat_stub.mod $(COMMON)/multicom_serial_string.mod

$(CVOBJ) cv_cvrms.o : cv_types.o cv_common.o
cv_posi_com.o cv_quaternion.o: cv_frames.o 
cv_proj.o : cv_drmsd.o
cv_qcomp.o : cv_quaternion.o
cv_cvrms.o : $(CVOBJ)
smcv.o : smcv_master.o cv_frames.o cv_quaternion.o cv_common.o
smcv_add.o : $(CVOBJ) cv_quaternion.o cv_types.o cv_common.o cv_frames.o
smcv_master.o : cv_common.o $(CVOBJ)
