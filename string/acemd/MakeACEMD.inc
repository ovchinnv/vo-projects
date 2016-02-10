LIB=libsmcv_acemd.so
LIBS+=$(LIB)
LINKOBJS=constants.o output.o parser.o cmd.o prm.o ivector.o ivector_list.o parselist.o files.o lu.o multidiag.o isort.o bestfit.o\
         sm_var.o sm_config.o ftsm_var.o smcv.o smcv_add.o smcv_master.o sm_util.o smcv_stats.o\
         cv_types.o cv_common.o cv_frames.o cv_quaternion.o \
         cv_angle_com.o cv_anglvec.o cv_posi_com.o cv_dihe_com.o cv_dist_com.o cv_drmsd.o cv_rmsd.o cv_proj.o cv_cvrms.o\
         cv_qcomp.o \
         plugin-smcv.o minimal_system.o smcv_acemd.o splines.o

default : acemd
acemd : all $(LIB)
$(LIB) : acemd.stamp MakeACEMD.inc
