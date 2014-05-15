# THIS MAKEFILE IS GENERATED AUTOMATICALLY FROM A MASTER SOURCE CODE TREE
ROOT=..
COMMON=$(ROOT)/scommon
CVOBJ=cv_posi_com.o cv_dist_com.o cv_angle_com.o cv_anglvec.o cv_dihe_com.o cv_rmsd.o cv_drmsd.o cv_proj.o cv_quaternion.o cv_qcomp.o
OBJS=cv_types.o cv_common.o cv_frames.o cv_cvrms.o smcv_master.o smcv.o smcv_add.o $(CVOBJ)
LOCALMACROS=-include source.defs -include $(COMMON)/source.defs 
LOCALDEPS=source.defs $(ROOT)/mpitype.def $(COMMON)/source.defs 
CHARMM_SOURCES=$(OBJS:%.o=%.src)
include $(ROOT)/Makefile.inc
#### RULES:####
charmm : $(CHARMM_SOURCES)


clean ::
	rm -f $(CHARMM_SOURCES) *tmp
