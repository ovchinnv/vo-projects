LIB=libstruna_namd.so
LIBS+=$(LIB)
LINKOBJS=constants.o output.o parser.o cmd.o prm.o\
         ivector.o ivector_list.o parselist.o files.o lu.o multidiag.o isort.o bestfit.o\
         sm_var.o sm_config.o smcv.o smcv_add.o smcv_master.o sm_util.o smcv_stats.o smcv_check_set.o\
         cv_types.o cv_common.o cv_frames.o cv_quaternion.o\
         cv_angle_com.o cv_anglvec.o cv_posi_com.o cv_dihe_com.o cv_dist_com.o cv_drmsd.o cv_rmsd.o cv_proj.o cv_cvrms.o\
         cv_qcomp.o system.o sm_namd.o splines.o rng.o clcg.o mpi_compat_stub.o psf.o ch_param.o mol_formats.o multicom_aux.o\
         psfatom.o state.o fieldbcdef.o formats.o atompar.o bondpar.o anglpar.o dihepar.o tlist.o angles.o bonds.o dihedrals.o\
         ftsm.o ftsm_util.o ftsm_io.o ftsm_voronoi.o ftsm_rex.o ftsm_rep.o ftsm_min.o ftsm_connect.o ftsm_compute.o ftsm_addatoms.o\
         ftsm_interpolate.o ftsm_stats.o ftsm_var.o ftsmv2_compute.o tsp.o\
         corio.o pdbio.o charmmio.o freeio.o chestio.o plot3dio.o\
         rvector.o rvector_list.o rsort.o util.o stats.o

default : namd
namd : $(LIB)
$(LIB) : $(LOBJS) MakeNAMD.inc

GlobalMasterSM.h:GlobalMasterSM.H $(LOCALDEPS) $(GLOBALDEPS)
	$(FPP) $(FPPFLAGS) $(DEFS) $(LOCALDEFS) $(GLOBALMACROS) $(LOCALMACROS) -P $< > $@
	cat source.defs >> GlobalMasterSM.h


