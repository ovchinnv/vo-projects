LIB=libstruna.so libstruna.a
LIBS+=$(LIB)
LINKOBJS=\
 atompar.o \
 bondpar.o \
 anglpar.o \
 dihepar.o \
 tlist.o \
 bonds.o \
 angles.o \
 dihedrals.o \
 pdbio.o \
 charmmio.o \
 chestio.o \
 corio.o \
 freeio.o \
 plot3dio.o \
 constants.o \
 output.o \
 parser.o \
 cmd.o prm.o \
 ivector.o \
 ivector_list.o \
 isort.o \
 rvector.o \
 rvector_list.o \
 rsort.o \
 parselist.o \
 files.o \
 lu.o \
 multidiag.o \
 bestfit.o \
 sm_var.o \
 sm_config.o \
 smcv.o \
 smcv_add.o \
 smcv_master.o \
 sm_util.o \
 smcv_stats.o \
 cv_types.o \
 cv_common.o \
 cv_frames.o \
 cv_quaternion.o \
 cv_qcomp.o \
 cv_angle_com.o \
 cv_anglvec.o \
 cv_posi_com.o \
 cv_dihe_com.o \
 cv_dist_com.o \
 cv_drmsd.o \
 cv_rmsd.o \
 cv_proj.o \
 cv_cvrms.o \
 ftsm_var.o ftsm.o \
 ftsm_addatoms.o \
 ftsm_compute.o \
 ftsm_connect.o \
 tsp.o \
 ftsm_interpolate.o \
 ftsm_io.o \
 ftsm_min.o \
 ftsm_rep.o \
 ftsm_rex.o \
 ftsm_stats.o \
 ftsm_util.o \
 ftsmv2_compute.o \
 ftsm_voronoi.o \
 minimal_system.o \
 splines.o \
 rng.o clcg.o \
 mpi_compat_stub.o \
 system.o \
 psf.o \
 ch_param.o \
 mol_formats.o \
 smcv_check_set.o \
 psfatom.o \
 state.o \
 fieldbcdef.o \
 formats.o \
 stats.o \
 util.o


default : plugin
plugin : $(LIB)
$(LIB) : $(LOBJS) Makefile.inc

