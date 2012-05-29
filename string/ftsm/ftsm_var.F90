/*COORDINATES AND MASSES:*/
/*
#ifdef __IMPNONE
#undef __IMPNONE
#endif
#define __IMPNONE
*/
! **********************************************************************!
! This source file was was generated automatically from a master source !
! code tree, which may not be distributed with this code if the !
! distributor has a proprietary compilation procedure (e.g. CHARMM) !
! If you edit this file (rather than the master source file) !
! your changes will be lost if another pull from the master tree occurs.!
! In case you are wondering why, this approach makes it possible for !
! me to have the same master source code interfaced with different !
! applications (some of which are written in a way that is quite far !
! from being object-oriented) at the source level. !
! **********************************************************************!
!
! FTSM_VAR.MOD
!
! VARIABLES FOR THE FINITE TEMPERATURE STRING METHOD
!**CHARMM_ONLY**!##IF STRINGM
!
      module ftsm_var
!
! note that SMCV and FTSM should not be used simultaneously
! see module smcv_common for the description of the variables below
!
      use sm_var, only: &
     & nstring, mestring, &
     & repa_initialized, smcv_initialized, &
     & linear, spline, bspline, interp_method, iterations, def, &
     & dst_cutoff, &
     & ds, curv, & ! arclength and curvature
     & stat_iteration_counter, &
     & output_rmsd0, &
     & output_arclength, &
     & output_curvature, &
     & output_fe, &
     & output_forces, &
     & output_rex_map, &
     & output_rex_log, &
!
     & stat_initialized, &
!
     & forces_fname, rex_fname, &
     & rmsd0_fname, s_fname, fe_fname, c_fname, &
!
     & rmsd0_funit, s_funit, &
     & fe_funit, c_funit, &
     & forces_funit, rex_funit, &
!
     & rform, sform, feform, cform, fform, rxlform, &
!
     & forces_flen, rmsd0_flen, s_flen, fe_flen, c_flen, rex_flen
!
      use sm_config, only: &
     & ftsm_on, repa_on, restrained_on, &
     & stat_on, &
     & string_noprint, restraint_force_on, &
     & calc_bestfit_grad_para, &
     & evolve_freq, stat_freq, &
     & restrained_eq_steps, restrained_eq0, evolve_nskip, &
     & olditeration, &
     & evolve_aver_on, & ! in this context, whether the evolution corresponds to averaging the simulation structure
! ! when false, evolve_expo_on is used, which implies a fixed exponential filter width
     & evolve_expo_on, &
     & evolve_expo_mem, &
     & finite_difference_d, &
     & parallel_tolerance, &
     & repl_x_on, &
     & repl_x_freq, &
     & rextime_offset
!
       real*8, allocatable, save :: fe(:), feav(:) ! free energy arrays
       real*8, save :: avforce(2) ! average parallel (1) and (2) perpendicular force;
       real*8, save :: repl_x_temp ! temperature for replica exchange
       integer, save :: num_evolve_samples=0 ! number of samples in the averaged image
       integer, save :: max_evolve_samples=0 ! maximum number of allowed samples before wraparound (if > 0)
       integer, save :: num_fe_samples=0 ! number of f.e. curves in average
       integer, save :: num_force_samples=0 ! number of force samples
       logical, save :: ftsm_initialized=.false.
       logical, save :: evolve_ftsm_on=.false., & ! is string evolution on?
     & update_on=.false. ! image updating on?
       logical, save :: output_centers=.false. ! center of the transition tube (these are coordinate files)
       logical, save :: qdiffrot=.false. ! whether the orientation atoms are different from forcing atoms, in which case store them separately
       logical, save :: qorient=.false. ! whether images are to be oriented in the bestfit sense
       logical, save :: proj_on=.false. ! whether the string evolves only along the direction perpendicular to itself
       character(len=80), save :: centers_fname=''
       integer, save :: centers_funit=-1, centers_flen=0
       character(len=80), save :: cenform
!
       integer, save :: norient=0, nforced=0, & ! orientation and forcing atoms
     & nany=0, nboth=0 ! any atoms, overlapping atoms
       integer, save :: update_freq=0 ! frequency for performing image update
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! coordinate variables below
       real*8, save, pointer :: &
     & r_f(:,:,:), &! forcing
     & r_o(:,:,:) ! orientation
       real*8, save, pointer :: rcom(:,:)
! weights
       real*8, save, pointer :: orientWeights(:), forcedWeights(:)
! indices
       integer, save, pointer :: iatom_o(:), iatom_f(:), iatom_a(:)
       integer, save, pointer :: iatom_both(:,:) ! index pairs that correspond to the same PSF atoms
! define indices into r array:
       integer, parameter, public :: left=1, & !
     & center=2, & !
     & right=3, & !
     & left_old=4, & ! old structures are saved to provide smooth transitions after evolution
     & center_old=5, & !
     & right_old=6, & !
     & center_new=7, & !
     & ref=8, & !
     & rave=9, & !
     & dummy=10, & !
     & instant=11, & !
     & vpar=12, & !
     & vperp=13, & !
     & fpar=14, & !
     & fperp=15, & !
     & left_rot=16, & !
     & center_rot=17, & !
     & right_rot=18, & !
     & left_cur=19, &
     & center_cur=20, &
     & right_cur=21
!
       real*8, save :: krms=0d0, kpara=0d0, kperp=0d0, &
     & dpar0=0.5d0,dperp0=0d0,drms0=0d0,dpar,dperp,drms
!
       integer, parameter, public :: num_sets=21 ! num of parameters above
       integer, save :: MPI_RTMD_TYPE, MPI_RTMD_TYPE_
!
       character(len=8), parameter, public :: real_format='(E23.15)'
       character(len=5), parameter, public :: int_format='(I10)'
!
      end module ftsm_var
!
!**CHARMM_ONLY**!##ENDIF
