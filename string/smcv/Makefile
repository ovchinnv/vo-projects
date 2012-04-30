ROOT=../..
CVOBJ=cv_posi_com.o cv_dist_com.o cv_angle_com.o cv_anglvec.o cv_dihe_com.o cv_rmsd.o cv_drmsd.o cv_proj.o cv_quaternion.o cv_qcomp.o
OBJS=cv_types.o cv_common.o cv_frames.o cv_cvrms.o $(CVOBJ)

COMMON=..
INCLUDES= -I $(COMMON) -I $(ROOT)/modules

LOCALMACROS=-include source.defs -include ../source.defs
LOCALDEPS=source.defs ../source.defs

include $(ROOT)/Makefile.inc

cv_common.o cv_frames.o: sm_common_vars output vectors

$(CVOBJ) cv_cvrms.o : cv_types.o cv_common.o multicom_aux psf

cv_posi_com.o cv_quaternion.o: cv_frames.o 
cv_frames.o : bestfit multicom_aux multicom_aux
cv_common.o cv_dihe_com.o cv_cvrms.o: constants
cv_rmsd.o cv_drmsd.o cv_proj.o: bestfit
cv_proj.o : cv_drmsd.o
cv_qcomp.o : cv_quaternion.o
cv_cvrms.o : $(CVOBJ)

sm_common_vars:
	$(MAKE) -C $(COMMON) sm_var.o sm_config.o 

charmm : $(CHARMM_SOURCES)

.PHONY : sm_common_vars
