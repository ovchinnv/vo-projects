/*#define __WRN(__WHO,__MSG) write(0,*) 'WARNING FROM: ',__WHO,': ',__MSG*/
/*#define __PRINT(__MSG) write(0,'(A)') __MSG*/
/*COORDINATES AND MASSES:*/
/*#define __INDX(__STR, __STRLEN, __TEST, __TESTLEN)  index(__STR(1:min(__STRLEN,len(__STR))),__TEST(1:min(__TESTLEN,len(__TEST))))*/
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
!CHARMM Element source/stringm/sm0k.src $Revision: 1.5 $
!
! zero-temperature string code
! documentation will be provided in stringm.doc
!
!**CHARMM_ONLY**!##IF STRINGM
!
      module sm0k ! string-method-at-0-K
      private
!
! VARIABLES
!ccccc initialization flag
      logical, public, save :: sm0k_initialized=.false.
!ccccc number of replicas on the string
      integer, save :: nstring=-1
      integer, save :: mestring=-1
!ccccc GENERAL VARS
      logical, save :: repa_initialized
      integer, save :: norient, nmove
! interpolation methods
      integer, parameter :: linear=1, spline=2, bspline=3, dst=4, &
     & linear_exact=5
!
      integer, save :: interp_method=0,orient=0
      logical, save :: qstat_orient=.false.
      integer, save :: orient_mass=0, repa_mass=0
      integer, save :: iterations=1 ! maximum interpolation iterations
      real*8, save :: def=1.1d0 ! interpolation tolerance
      real*8, save :: dst_cutoff=1.0d0 ! wavenumber truncation parameter for DST
!
! arrays
      real*8, save, allocatable :: &
     & rcurrent_m(:,:), rref_o(:,:), rcurrent_o(:,:)
! arclength and curvature
      real*8, save, allocatable :: ds(:), curv(:) ! unavailable at first iteration
! orientation weights -- for 0-temp. string
      real*8, save, allocatable :: orientWeights(:), repaWeights(:)
      integer, save, allocatable :: iatom_o(:), iatom_m(:), &
     & iatom_free_o(:), iatom_free_m(:)
      integer, save, pointer :: iatom_f(:)
      logical, save, allocatable :: fixed_o(:), fixed_m(:), fixed_s(:)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccc STATISTICS VARS cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      integer, save :: stat_iteration_counter=0 ! how many times 'stat' has been called
      logical, save :: output_rmsd0=.false., &
     & output_dsdt=.false., &
     & output_arclength=.false., & ! output options
     & output_curvature=.false., &
     & output_rmsd_ave=.false. ! rmsd wrt to the average structure
      logical, save :: stat_rmsd_mass=.false. ! should the statistics routine do mass-weighted RMSD?
      logical, save :: stat_initialized=.false.
!
      character(len=80), save :: rmsd0_fname='', dsdt_fname='', s_fname='', &
     & rmsd_ave_fname='', c_fname='' ! output names
!
      integer, save :: nstat
      integer, save :: rmsd0_funit=-1, dsdt_funit=-1, s_funit=-1, &
     & rmsd_ave_funit=-1, c_funit=-1
      integer, save :: num_average_samples=0 ! number of samples in the average set
      integer, save :: rmsd0_flen=0, dsdt_flen=0, s_flen=0, rmsd_ave_flen=0, c_flen
! arrays
      real*8, save, allocatable, dimension(:,:) :: &
     & rold_s, rave_s, rcurrent_s, rcomp_s, rold_o, rave_o, rcomp_o
      real*8, save, allocatable :: statWeights(:)
      integer, save, allocatable :: iatom_s(:), iatom_free_s(:)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! SUBROUTINES
!
      public sm0k_main
      public sm0k_repa
      public sm0k_stat
!
      contains
!
      SUBROUTINE sm0k_main(COMLYN,COMLEN)
!----------------------------------------------------------------------
! command parser for the 0K string
!----------------------------------------------------------------------
      use sm_config, only : stat_on, stat_freq, repa_on, repa_freq ! for communication with MINI
!
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use psf
      use system, only : r, rcomp, m, bfactor, occupancy
      use parser
!
      implicit none
!
      CHARACTER(LEN=*) :: COMLYN
      INTEGER :: COMLEN
!
! local variables
      character(len=8) :: keyword
      character(len=11) :: whoami
!
  character(len=80) :: msg___
!
! declare functions here
!
      data whoami /' SM0K_MAIN>'/
!
      keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (( keyword(1:4).eq.'INIT'(1:4) )) then
        call sm0k_init()
        return
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (.not.sm0k_initialized) then
        call sm0k_init()
      endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (( keyword(1:4).eq.'INTE'(1:4) )) then
        call sm0k_interpolate(comlyn, comlen)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'REPA'(1:4) )) then
       if (comlen.gt.0) then ! this is an initialization call!
        call sm0k_repa_init(comlyn, comlen)
       else
        call sm0k_repa(0)! repa routine will reparametrize main coords
       endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'STAT'(1:4) )) then
       if (comlen.gt.0) then ! this is an initialization call!
        call sm0k_stat_init(comlyn, comlen)
       else
        call sm0k_stat(0) ! compute statistics from main coordinates
       endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'MINI'(1:4) )) then
! string minimization
! SD routine will be called;
! SD is the only minimizer allowed, other minimizers removed below
! setting "repa_on" to true so that
! SD knows to call reparametrization
! other options are let through -- use at your risk!
! delete ABNR, POWE, CONJ, CGSD
!cccccccccccccccccc reparametrization option cccccccccccccccccccccc
       repa_freq=atoi(get_remove_parameter(comlyn, 'REPF', comlen), -1)
       if (repa_freq.le.0) then
        repa_on=.false.
        WRITE (msg___,'(/,2A,/,2A,/,2A/)') &
     & whoami,' STRING METHOD ENABLED, BUT', &
     & whoami,' REPARAMETRIZATION FREQUENCY ZERO OR UNSPECIFIED.', &
     & whoami,' REPARAMETRIZATION WILL NOT BE DONE.' ; call plainmessage(msg___,3)
       else
        repa_on=.true.
        WRITE (msg___,'(/,2A,/,2A,I7,A/)') &
     & whoami,' STRING METHOD ENABLED.', &
     & whoami,' WILL REPARAMETRIZE AFTER EVERY ', &
     & repa_freq,' MINIMIZATION ITERATIONS' ; call plainmessage(msg___,3)
       endif ! repa_freq
!cccccccccccccccccc statistics output option cccccccccccccccccccccc
       if (repa_on) then ! boolly, it makes sense to output string statistics only when reparametrization is enabled
! if you want to follow the unparametrized dynamics, just set maxiter to 0 in the repa setup call
        stat_freq=atoi(get_remove_parameter(comlyn, 'STAF', comlen), -1)
        if (stat_freq.le.0) then
        stat_on=.false.
        WRITE (msg___,'(/,2A,/,2A/)') &
     & whoami,' STATISTICS OUTPUT FREQUENCY NOT SPECIFIED.', &
     & whoami,' STATISTICS WILL NOT BE OUTPUT.' ; call plainmessage(msg___,3)
        else
         stat_on=.true.
         WRITE (msg___,'(/,2A,I6,A/)') &
     & whoami,' WILL OUTPUT STATISTICS AFTER EVERY ', &
     & stat_freq,' REPARAMETRIZATION ITERATIONS' ; call plainmessage(msg___,3)
        stat_freq=stat_freq*repa_freq
        endif ! stat_freq
       else ! repa_on
        WRITE (msg___,'(/,2A,/,2A/)') &
     & whoami,' STATISTICS OUTPUT REQUIRES REPARAMETRIZATION', &
     & whoami,' (DISABLED). STATISTICS WILL NOT BE OUTPUT.' ; call plainmessage(msg___,3)
        stat_on=.false.
       endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      else
            write(msg___,*)'UNRECOGNIZED SUBCOMMAND: ',keyword;call warning(whoami, msg___, 0)
      endif
!
      end subroutine sm0k_main
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sm0k_init()
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use parser
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use mpi
!
      implicit none
 character(len=80) :: msg___
!
      logical :: qroot, qslave
      integer :: ierror
      character(len=11) :: whoami
!
      data whoami /' SM0K_INIT>'/
!
! do a basic communicator check:
      if (ME_LOCAL.eq.0.and.ME_STRNG.eq.MPI_UNDEFINED) then
        write(msg___, 111) whoami, ME_GLOBAL, whoami ; call plainmessage(msg___)
 111 FORMAT(A, ' WORLD REPLICA ',I5, ' HAS ZERO GROUP ID', &
     & /,A,' BUT INVALID STRING ID (MAY BE OK).')
      elseif (ME_STRNG.ne.MPI_UNDEFINED.and. &
     & (ME_LOCAL.ne.0.or.MPI_COMM_LOCAL.eq.MPI_COMM_NULL)) then
        write(msg___, 111) whoami, ME_GLOBAL, whoami ; call plainmessage(msg___)
 112 FORMAT(A, ' WORLD REPLICA ',I5, ' HAS A VALID STRING ID', &
     & /,A,' BUT A NONZERO GROUP ID. ABORTING.')
       return
      endif
!
      qroot=ME_STRNG.ne.MPI_UNDEFINED
      qslave=ME_LOCAL.ne.MPI_UNDEFINED ! (also includes roots)
!
      if (sm0k_initialized) then
       if (qroot) then
        if (ME_STRNG.eq.0) then
          write(msg___,'(2A)') &
     & whoami, ' SM0K ALREADY INITIALIZED. CALL "DONE" TO CLEAN UP.' ; call plainmessage(msg___)
        endif
       endif
       return
      endif
!
      nstring=1 ! safe (hopefully) default
      mestring=-1 ! safe (hopefully) default
!
      if (qroot) then
        nstring=SIZE_STRNG
        mestring=ME_STRNG
      endif
      call mpi_bcast(nstring,1,MPI_INTEGER,0,MPI_COMM_LOCAL,ierror)
      call mpi_bcast(mestring,1,MPI_INTEGER,0,MPI_COMM_LOCAL,ierror)
!
      if (qroot) then
        if (ME_STRNG.eq.0) then
          write(msg___,'(2A,I5, A)') &
     & whoami, ' FOUND ',nstring,' REPLICAS.' ; call plainmessage(msg___)
        endif
      endif
! store fixed atom indices
      allocate(iatom_f(0))
!
      sm0k_initialized=.true.
!
      end subroutine sm0k_init
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sm0k_done()
      use parser
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use mpi
      implicit none
 character(len=80) :: msg___
!
      character(len=11) :: whoami
!
      data whoami /' SM0K_DONE>'/
      if (MPI_COMM_STRNG.ne.MPI_COMM_NULL.and.ME_STRNG.eq.0) then
       write(msg___,'(2A,I5, A)') whoami, ' CLEANING UP.' ; call plainmessage(msg___)
      endif
      sm0k_initialized=.false.
      nstring=-1
      mestring=-1
!
      if (associated(iatom_f)) deallocate(iatom_f)
!
      if (allocated(ds)) deallocate(ds)
      if (allocated(curv)) deallocate(curv)
!
      end subroutine sm0k_done
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sm0k_repa_init(COMLYN, COMLEN)
! initialize string reparametrization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use psf
      use parser
      use system, only : r, rcomp, m, bfactor, occupancy
      use system, only : system_getind
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
 character(len=80) :: msg___
!
      CHARACTER(LEN=*) :: COMLYN
      INTEGER :: COMLEN
!
      character(len=22) :: methods(5)
      character(len=16) :: whoami
      character(len=8) :: keyword
      data methods &
     & / 'LINEAR','CUBIC SPLINE','B-SPLINE','DST','LINEAR EXACT'/
! selection array
      integer :: i ,j, mlen
      integer :: qlinear, qspline, qbspline, qdst, qlinear_exact
!
  integer :: isele, i__, iend; integer, pointer::iselct(:);character(LEN=20)::word__
      integer, pointer :: kselct(:) ! for swapping iselct/jselct
!
      integer, pointer :: jselct(:)
      integer :: nslct
!
      data whoami /' SM0K_REPA_INIT>'/
!
      if (.not.sm0k_initialized) call sm0k_init()
!
! reset variables
      qspline=0
      qbspline=0
      qlinear=0
      qdst=0
      qlinear_exact=0
      dst_cutoff=0.0
      interp_method=0
      orient=0
      orient_mass=0 ! will the orientation use mass weighting?
      repa_mass=0 ! will the reparametrization use mass weighting
      nmove=0
      norient=0
      num_average_samples=0
!
! deallocate arrays
      if (allocated(fixed_o)) deallocate(fixed_o) ! flags
      if (allocated(fixed_m)) deallocate(fixed_m) ! flags
      if (allocated(iatom_o)) deallocate(iatom_o)
      if (allocated(iatom_m)) deallocate(iatom_m)
      if (allocated(iatom_free_o)) deallocate(iatom_free_o)
      if (allocated(iatom_free_m)) deallocate(iatom_free_m)
      if (allocated(rref_o)) deallocate(rref_o)
      if (allocated(rcurrent_o)) deallocate(rcurrent_o)
      if (allocated(rcurrent_m)) deallocate(rcurrent_m)
      if (allocated(orientWeights)) deallocate(orientWeights)
      if (allocated(repaWeights)) deallocate(repaWeights)
!
! deallocate arclength+curavture arrays
      if (allocated(ds)) deallocate(ds)
! initialize curvature array
      if (allocated(curv)) deallocate(curv)
!
      repa_initialized=.false.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (remove_tag(comlyn,'LINE',comlen).gt.0) then
       qlinear=1
       interp_method=linear
      endif
      if ((remove_tag(comlyn,'CSPL',comlen).gt.0).or. &
     & (remove_tag(comlyn,'SPLI',comlen).gt.0)) then
       qspline=1
       interp_method=spline
      endif
      if (remove_tag(comlyn,'BSPL',comlen).gt.0) then
       qbspline=1
       interp_method=bspline
      endif
      if (remove_tag(comlyn,'DST',comlen).gt.0) then
       qdst=1
       interp_method=dst
! did the user specify filter cutoff?
       dst_cutoff=atof(get_remove_parameter(comlyn, 'WNCT', comlen), -1.0d0)
       if (dst_cutoff.lt.0.0) then
        write(msg___,664) whoami, whoami ; call plainmessage(msg___,3)
 664 FORMAT(A,' DST REQUESTED BUT FILTER CUTOFF', &
     & A, ' NOT SPECIFIED.',/,' WILL USE 0.500')
        dst_cutoff=0.5
       endif
      endif
      if (remove_tag(comlyn,'LIN2',comlen).gt.0) then
       qlinear_exact=1
       interp_method=linear_exact
      endif
!ccccccc CHECK FOR MULTIPLE OPTIONS
      if ((qspline+qlinear+qbspline+qdst+qlinear_exact) .eq. 0) then
       write(msg___,665) whoami, whoami ; call plainmessage(msg___,3)
 665 FORMAT(A,' INTERPOLATION METHOD NOT SPECIFIED.',/, &
     & A,' WILL USE LINEAR INTERPOLATION.')
       interp_method=linear
      elseif ((qdst+qspline+qlinear+qbspline+qlinear_exact) .gt. 1)then
       call warning(whoami, 'TOO MANY INTERPOLATION OPTIONS.', 0)
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! did the user specify an interpolation tolerance?
      if (interp_method.ne.linear_exact) then ! options below are invalid for exact interpolation
       def=atof(get_remove_parameter(comlyn, 'DEFI', comlen), 1.1d0)
       if (def.lt.1.0) then
         call warning(whoami, 'INTERPOLATION TOLERANCE MUST BE >= 1. EXITING.', 0)
         return
       endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! did the user specify a maximum number of iterations?
       iterations=atoi(get_remove_parameter(comlyn, 'ITER', comlen), 10)
      else
       def=0d0
       iterations=0
      endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! check for orientation options: must appear after 'ORIE'
! process selection
!
!
      i=remove_tag(comlyn,'ORIE',comlen) ! find the position of `intc'
      if (i.gt.0) then ! only if the ORIE directive exists
! selection text taken from corman.src
       orient=1
       j=find_tag(comlyn, 'SELE', comlen)
       if (j.gt.0.and.j.lt.i) then ! sele occurs before orie
        call warning(whoami, 'ATOM SELECTION MUST BE SPECIFIED AFTER ORIE.', 0)
        return
       endif
!
          nullify(iselct)
          isele=find_tag(comlyn, 'SELE', comlen)
          msg___=comlyn(isele:comlen) ! part of string that begins with the selection
          i__=comlen-isele+1
          iend=find_tag(msg___, 'END', i__) ! location of selection termination
          iend=iend-1+isele ! index into comlyn starting from 1
          msg___=comlyn(isele:iend) ! part of string that begins with the selection and ends before ' END'
          word__=pop_string(msg___) ! remove first word (which we know is 'SELE*') from msg___
! process selection:
          nullify(iselct)
          iselct=>system_getind(msg___)
! remove selection string from command line:
          msg___=comlyn(iend:comlen) ! command line starting with 'END' (see above)
          word__=pop_string(msg___) ! remove 'END*' e.g. 'ENDING' is ok too
          comlyn(isele:isele)=' ';
          comlyn(isele+1:)=msg___ ! selection has been removed from command
          comlen=len_trim(comlyn)
      endif ! orie specified
!
! check for mass weighting
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      j=1
      do while (j.gt.0)
       mlen=comlen
       j=remove_tag(comlyn,'MASS',comlen) ! s
       if ( (orient.eq.1).and.(j.ge.i) ) then
        orient_mass=1
       else if (j.gt.0) then
        repa_mass=1
        i=i-(mlen-comlen) ! how much the string has shrunk
       endif
      enddo
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! check which atoms will be moved by reparametrization (default: all atoms)
! process selection
      i=remove_tag(comlyn,'MOVE',comlen) !
      if (i.gt.0) then ! only if the MOVE directive exists
       j=find_tag(comlyn, 'SELE', comlen)
       if (j.gt.0.and.j.lt.i) then ! sele occurs before orie
        call warning(whoami, 'ATOM SELECTION MUST BE SPECIFIED AFTER MOVE.', 0)
        return
       endif
!
       kselct=>iselct ; nullify(iselct)
          nullify(iselct)
          isele=find_tag(comlyn, 'SELE', comlen)
          msg___=comlyn(isele:comlen) ! part of string that begins with the selection
          i__=comlen-isele+1
          iend=find_tag(msg___, 'END', i__) ! location of selection termination
          iend=iend-1+isele ! index into comlyn starting from 1
          msg___=comlyn(isele:iend) ! part of string that begins with the selection and ends before ' END'
          word__=pop_string(msg___) ! remove first word (which we know is 'SELE*') from msg___
! process selection:
          nullify(iselct)
          iselct=>system_getind(msg___)
! remove selection string from command line:
          msg___=comlyn(iend:comlen) ! command line starting with 'END' (see above)
          word__=pop_string(msg___) ! remove 'END*' e.g. 'ENDING' is ok too
          comlyn(isele:isele)=' ';
          comlyn(isele+1:)=msg___ ! selection has been removed from command
          comlen=len_trim(comlyn)
       jselct=>iselct
       iselct=>kselct
      else
       jselct=>system_getind('ALL') ! select all atoms by default
      endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! now process the selections
! 3/11 provide compatibiliy with fixed atoms and nominal compatibility with TSM backlists
!
      if (associated(iselct)) then ; norient=size(iselct) ; else ; norient=0 ; endif
!
      if (norient.eq.0) then
        call warning(whoami, 'NO ATOMS SELECTED FOR ORIENTATION. WILL NOT ORIENT.', 0)
        orient=0 ! invalid or missing selection for orientation
       elseif (norient.lt.3) then
        call warning(whoami, 'FEWER THAN FOUR ATOMS SELECTED FOR ORIENTATION. WILL NOT ORIENT.', 0)
        orient=0
      endif
!
      if (associated(jselct)) then ; nmove=size(jselct) ; else ; nmove=0 ; endif
!
! allocate space for various arrays atom array
      allocate(iatom_o(norient), iatom_m(nmove), &
     & iatom_free_o(norient), iatom_free_m(nmove), & ! these index lists are for use with minimizer arrays (fixed atoms removed)
     & fixed_o(norient), fixed_m(nmove), & ! these flags indicate that the atom will be absent from minimizer array
     & rref_o(norient,3),rcurrent_o(norient,3), &
     & rcurrent_m(nmove,3), &
     & orientWeights(norient), repaWeights(nmove))
!
      iatom_o=0; iatom_m=0;
      rref_o=0d0
      rcurrent_o=0d0
      rcurrent_m=0d0
      orientWeights=1d0
      repaWeights=1d0
! build various index arrays
!
! currently, no support for fixed atoms -- assume their number is zero
! orientation
      if (associated(iselct)) then
       iatom_o=iselct
       iatom_free_o=iatom_o*3-2 ! gives 1,4,7...
       fixed_o=.false.
       deallocate(iselct)
      endif
! moving
      if (associated(jselct)) then
       iatom_m=jselct
       iatom_free_m=iatom_m*3-2
       fixed_m=.false.
       deallocate(jselct)
      endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (nmove.eq.0) then ! this will occur if all of the atoms selected are fixed
        call warning(whoami, 'NO ATOMS SELECTED FOR REPARAMETRIZATION. CANNOT CONTINUE.', 0)
        if (associated(iselct)) deallocate(iselct)
        if (associated(jselct)) deallocate(jselct)
        return
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! done with selections
!!!!!!!!! process mass-weighting
      if (orient_mass.eq.1) then
       do i=1,norient
        orientWeights(i)=m(iatom_o(i))*orientWeights(i) ! these weights are for the Best-Fit routine
! orientWeights(i)=sqrt(amass(iatom_o(i)))*orientWeights(i)
       enddo
!
       do i=1, nmove
! repaWeights(i)=sqrt(amass(iatom_m(i)))*repaWeights(i) ! these weights are essentially for multiplications
        repaWeights(i)=m(iatom_m(i))*repaWeights(i)
       enddo
      endif
!!!!!!
! normalize orientWeights;
      orientWeights=orientWeights/sum(orientWeights)
! unnecessary -- repaWeights are normalized in interpolation routine
      repaWeights=repaWeights/sum(repaWeights)
!
! print summary
!!!!!! reparametrization summary !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      mlen=len(methods(interp_method))
mlen=min(max(0,mlen),len(methods(interp_method)));methods(interp_method)(mlen+1:)='';call adjustleft(methods(interp_method),(/' ',tab/));mlen=len_trim(methods(interp_method))
      if (interp_method.eq.linear_exact) then
        write(msg___,666) whoami,methods(interp_method)(1:mlen)
      else
        write(msg___,667) whoami,methods(interp_method)(1:mlen),whoami, &
     & def
      endif
      call plainmessage(msg___,3)
 666 format(A,' WILL REPARAMETRIZE STRING USING ',A,' INTERPOLATION')
 667 format(A,' WILL REPARAMETRIZE STRING USING ',A,/, &
     &A,' INTERPOLATION TO WITHIN MAX(DS)/MIN(DS) < ',F7.3,' TOLERANCE')
      if (iterations.gt.0) then ; write(msg___,668) whoami, iterations ; call plainmessage(msg___,3) ; endif
 668 format(A,' WITH A MAXIMUM OF ',I3,' ITERATIONS')
      if(interp_method.eq.dst) then ; write(msg___,6680) whoami,dst_cutoff*100.0 ; call plainmessage(msg___,3) ; endif
 6680 format(A,' DST INTERPOLATION WILL USE THE LOWER ',F8.4, &
     & '% OF WAVENUMBERS')
!
      write(keyword,'(I8)') nmove ; call plainmessage(msg___,3)
      mlen=len(keyword)
      mlen=min(max(0,mlen),len(keyword));keyword(mlen+1:)='';call adjustleft(keyword,(/' ',tab/));mlen=len_trim(keyword)
      write(msg___,669) whoami, keyword(1:mlen) ; call plainmessage(msg___,3)
 669 format(A,' INTERPOLATION IS BASED ON ',A,' CARTESIAN COORDINATES')
      if (repa_mass.eq.1) then ; write(msg___,672) whoami ; call plainmessage(msg___,3) ; endif
 672 format(A,' INTERPOLATION WILL USE MASS-WEIGHTING')
!
!!!!!! orientation summary !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (orient.eq.1) then
        write(keyword,'(I8)') norient ; call plainmessage(msg___,3)
        mlen=len(keyword)
        mlen=min(max(0,mlen),len(keyword));keyword(mlen+1:)='';call adjustleft(keyword,(/' ',tab/));mlen=len_trim(keyword)
        write(msg___,670) whoami, keyword(1:mlen) ; call plainmessage(msg___,3)
 670 format(A,' STRING WILL BE ORIENTED BASED ON ',A,' ATOMS')
        if (orient_mass.eq.1) then ; write(msg___,671) whoami ; call plainmessage(msg___,3) ; endif
 671 format(A,' ORIENTATION WILL USE MASS-WEIGHTING')
      endif ! orient
!
! initialize arclength array
      if (.not.allocated(ds)) then
       allocate(ds(nstring-1))
       ds=0.0d0
      endif
! initialize curvature array
      if (.not.allocated(curv)) then
       allocate(curv(nstring-2))
       curv=0.0d0
      endif
!
      repa_initialized=.true.
!
! deallocate temporary variables
      if (associated(iselct)) deallocate(iselct)
      if (associated(jselct)) deallocate(jselct)
!
      end subroutine sm0k_repa_init
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sm0k_repa(n,var)
! zero-temperature:reparametrization is based
! on Cartesian coordinates
! note: var is assumed to contain coordinates of free atoms (CHARMM SD minimizer convention) in [x1,y1,z1,x2,y2,z2...] triplet format
!
      use bestfit
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use constants
      use psf
      use system, only : r, rcomp, m, bfactor, occupancy
      use parser
      use mpi
!
      implicit none
!ccccccccccccccccccccccccccccccccccccccc
      integer :: n ! NOTE: arrays as optional arguments are quite a dangerous feature of F90
      real*8, optional :: var(*) ! ideally, n should give the dimension of n, but we live in an imperfect world
!
!
      integer :: stat(MPI_STATUS_SIZE)
      integer :: ierror
      integer :: me, ncpu
      integer :: i, j
      real*8 :: t
      logical :: qroot, qslave, qmanual
! aux arrays
      real*8 :: weights(nmove,3) ! for reparametrization
      real*8 :: u(3,3)=RESHAPE( (/1,0,0,0,1,0,0,0,1/), (/3,3/) ) ! rotation matrix
      real*8 :: rcurrent_com(3)=(/0d0,0d0,0d0/) ! COM vector
!
      character(len=11) :: whoami
! interface to reparametrization routine
! needed because of assumed shape array below
!
      interface
        subroutine interp_driver_sci(rin,rout,wgt,n, &
     & interp_method,tol,max_iterations,d_arclength, curvature, &
     & dst_cutoff, dr,r_bc_0, r_bc_1)
        use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
        implicit none
        integer :: n
        real*8 :: rin(n), rout(n), wgt(n)
        integer, intent(in) :: interp_method
        integer :: max_iterations
        real*8 :: tol, d_arclength(:), curvature(:)
        real*8, optional :: dst_cutoff
        real*8, optional :: dr(n) ,r_bc_0(n), r_bc_1(n)
        end subroutine interp_driver_sci
!
        subroutine interp_driver_sci_root(rin,rout,wgt,n, &
     & interp_method,tol,max_iterations,d_arclength, curvature, &
     & dst_cutoff, dr,r_bc_0, r_bc_1)
        use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
        implicit none
        integer :: n
        real*8 :: rin(n), rout(n), wgt(n)
        integer, intent(in) :: interp_method
        integer :: max_iterations
        real*8 :: tol, d_arclength(:), curvature(:)
        real*8, optional :: dst_cutoff
        real*8, optional :: dr(n) ,r_bc_0(n), r_bc_1(n)
        end subroutine interp_driver_sci_root
!
        subroutine interp_linear_exact(rin,rout,wgt,n, &
     & d_arclength, curvature, &
     & drout, &
     & r_bc_0, r_bc_1)
        integer :: n
        real*8 :: rin(n), rout(n), wgt(n)
        real*8 :: d_arclength(:), curvature(:)
        real*8, optional :: drout(n) ! optional computation of tangent
        real*8 , optional :: r_bc_0(n), r_bc_1(n) ! optional fixed bc data
       end subroutine interp_linear_exact
!
      end interface
!
      data whoami /' SM0K_REPA>'/
!
      if (.not.sm0k_initialized) call sm0k_init()
      qmanual=(n.eq.0)
!
      qroot =MPI_COMM_STRNG.ne.MPI_COMM_NULL.and.SIZE_STRNG.gt.1
      qslave=MPI_COMM_LOCAL.ne.MPI_COMM_NULL.and.SIZE_LOCAL.gt.1
!
! check if the user has made an initialization call
      if (.not.repa_initialized) then
       call warning(whoami, 'NO REPARAMETRIZATION OPTIONS SELECTED. NOTHING DONE.', 0)
       return
      endif
!
! load coordinates
! manual reparametrization: use PSF indices
      if (qmanual) then
! orientation
       do i=1,norient
        j=iatom_o(i)
        rcurrent_o(i,1)=r(1,j); rcurrent_o(i,2)=r(2,j); rcurrent_o(i,3)=r(3,j)
       enddo
! moving
       do i=1,nmove
        j=iatom_m(i)
        rcurrent_m(i,1)=r(1,j); rcurrent_m(i,2)=r(2,j); rcurrent_m(i,3)=r(3,j)
       enddo
      else ! `automatic', i.e. called from a minimization routine
! orientation atoms
       do i=1,norient
        if (fixed_o(i)) then ! grab coordinates from main coordinate array
         j=iatom_o(i)
         rcurrent_o(i,1)=r(1,j)
         rcurrent_o(i,2)=r(2,j)
         rcurrent_o(i,3)=r(3,j);
        else ! grab coordinates provided by minimizer
         j=iatom_free_o(i) ! x-index (y- z- indices follow)
         rcurrent_o(i,1)=var(j);j=j+1
         rcurrent_o(i,2)=var(j);j=j+1
         rcurrent_o(i,3)=var(j)
        endif
       enddo
! moving atoms
       do i=1,nmove
! here, we should be able to assume that the moving atom coords are passed in
        j=iatom_free_m(i) ! x-index (y- z- indices follow)
        rcurrent_m(i,1)=var(j);j=j+1
        rcurrent_m(i,2)=var(j);j=j+1
        rcurrent_m(i,3)=var(j)
       enddo
      endif ! qmanual
!
      if (orient.eq.1) then
! translate to centroid
       rcurrent_com=matmul(transpose(rcurrent_o), orientWeights)
!
       rcurrent_m(:,1)=rcurrent_m(:,1)-rcurrent_com(1)
       rcurrent_m(:,2)=rcurrent_m(:,2)-rcurrent_com(2)
       rcurrent_m(:,3)=rcurrent_m(:,3)-rcurrent_com(3)
!
       rcurrent_o(:,1)=rcurrent_o(:,1)-rcurrent_com(1)
       rcurrent_o(:,2)=rcurrent_o(:,2)-rcurrent_com(2)
       rcurrent_o(:,3)=rcurrent_o(:,3)-rcurrent_com(3)
      endif
!
! statistics: save current structure
      if (stat_initialized) then
        if (output_dsdt.or.output_rmsd_ave) then
! also need to deal with fixed atoms
         do i=1,nstat
          if (qmanual.or.fixed_s(i)) then ! grab coordinates from main coordinate array
           j=iatom_s(i)
           rold_s(i,1)=r(1,j)
           rold_s(i,2)=r(2,j)
           rold_s(i,3)=r(3,j)
          else
           j=iatom_free_s(i) ! x-index (y- z- indices follow)
           rold_s(i,1)=var(j);j=j+1
           rold_s(i,2)=var(j);j=j+1
           rold_s(i,3)=var(j)
          endif
         enddo
! stat orientation atoms
         if (qstat_orient) then
          rold_o=rcurrent_o ! COM-free
          rold_s(:,1)=rold_s(:,1) - rcurrent_com(1)
          rold_s(:,2)=rold_s(:,2) - rcurrent_com(2)
          rold_s(:,3)=rold_s(:,3) - rcurrent_com(3)
         endif
        endif ! dsdt or rmsd_ave
! update running average
        if (output_rmsd_ave) then
         t=1.0d0*num_average_samples/(num_average_samples+1)
         if (qstat_orient) then
          if (num_average_samples.gt.0) then
           call RMSBestFit(rold_o, rave_o, orientWeights, u)
           u=transpose(u)
           rold_o=matmul(rold_o, u)
           rold_s=matmul(rold_s, u)
          endif
          rave_o=t*rave_o+(1.0d0-t)*rold_o
         endif ! qstat_orient
         rave_s=t*rave_s+(1.0d0-t)*rold_s
         num_average_samples=num_average_samples+1
        endif
      endif
!
      if (qroot) then
       if (orient.eq.1) then
!ccccccccccc take care of orientation ccccccc
! send/receive orientation structure
        me=ME_STRNG
        ncpu=SIZE_STRNG
        if (me.gt.0) then
         call mpi_recv(rref_o,3*norient,MPI_REAL,me-1,1, &
     & MPI_COMM_STRNG, stat,ierror)
! orient rcurrent based on rref
! 12.09: using RTMD orientation routines
! no checking for undefined coordinates here
         call RMSBestFit(rcurrent_o,rref_o,orientWeights,u)
! transform current structure to overlap with reference
! (if orientation is off, u=I)
         u=transpose(u)
         rcurrent_o=matmul(rcurrent_o, u)
         rcurrent_m=matmul(rcurrent_m, u)
        endif
        if (me.lt.ncpu-1) then
         call mpi_send(rcurrent_o,norient*3,MPI_REAL,me+1,1, &
     & MPI_COMM_STRNG, ierror)
        endif ! me
       endif ! orient
!cccccccccccccc now call the approproate interpolation subroutine
       if (repa_mass.eq.1) then
        weights(:,1)=repaWeights(:)! need 3 sets for x-,y-,z- coords
        weights(:,2)=repaWeights(:)
        weights(:,3)=repaWeights(:)
       else
        weights=1.0d0
       endif
!
! call by name
! there is the following (compiler?) bug: the routine thinks that dr is present
! even though it is not -- can survive without this but this is disturbing
! call interp_driver_sci_root(RIN=rcurrent_m,ROUT=rcurrent_m,
       if (interp_method.eq.linear_exact) then
        call interp_linear_exact(RIN=rcurrent_m,ROUT=rcurrent_m, &
     & WGT=weights,N=3*nmove, D_ARCLENGTH=ds,CURVATURE=curv)
       else
        call interp_driver_sci(RIN=rcurrent_m,ROUT=rcurrent_m, &
     & WGT=weights,N=3*nmove, INTERP_METHOD=interp_method,TOL=def, &
     & MAX_ITERATIONS=iterations,D_ARCLENGTH=ds,CURVATURE=curv, &
     & DST_CUTOFF=dst_cutoff)
!
! call by argument order
! call interp_driver_sci(rcurrent_m,rcurrent_m,weights,3*nmove,
! & interp_method,def,iterations,ds,curv,dst_cutoff)
! call interp_driver_sci_root(rcurrent_m,rcurrent_m,weights,
! & 3*nmove,interp_method,def,iterations,ds,curv,dst_cutoff)
       endif
!
       if (orient.eq.1) then
        u=transpose(u)
        rcurrent_m=matmul(rcurrent_m, u) ! rotate back
! restore original COM
        rcurrent_m(:,1)=rcurrent_m(:,1)+rcurrent_com(1)
        rcurrent_m(:,2)=rcurrent_m(:,2)+rcurrent_com(2)
        rcurrent_m(:,3)=rcurrent_m(:,3)+rcurrent_com(3)
!
       endif ! orient
      endif ! root
!
! broadcast coordinates to slaves
      if (qslave) then
       call mpi_bcast(rcurrent_m,nmove*3,MPI_REAL,0,MPI_COMM_LOCAL,ierror)
      endif
!
! copy back moving atoms
!
      if (qmanual) then
       do i=1,nmove
        j=iatom_m(i)
        r(1,j)=rcurrent_m(i,1)
        r(2,j)=rcurrent_m(i,2)
        r(3,j)=rcurrent_m(i,3)
       enddo
      else
       do i=1,nmove
        if (fixed_m(i)) cycle ! do not return fixed atoms (should not be used)
        j=iatom_free_m(i) ! x-index (y- z- indices follow)
        var(j)=rcurrent_m(i,1);j=j+1
        var(j)=rcurrent_m(i,2);j=j+1
        var(j)=rcurrent_m(i,3)
       enddo
      endif
!
      end subroutine sm0k_repa
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sm0k_interpolate(comlyn, comlen)
! given a collection of n string replicas, interpolate onto a finer/coarser path
!
      use bestfit
!
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use mpi
!
      use parser
      use psf
      use constants
!
      use system, only : r, rcomp, m, bfactor, occupancy
!
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use system, only : system_getind
      use charmmio; use pdbio; use mol_formats
!
      implicit none
 character(len=80) :: msg___
!
      CHARACTER(LEN=*) :: COMLYN
      INTEGER :: COMLEN
!
! local declarations only: I do not want persistence, so no saving into modules
      integer, parameter :: linear=1, spline=2, bspline=3
      integer :: ifile, ofile, num_rep_in, num_rep_out
      integer :: len_cor_in, len_cor_out, &
     & length, mlen, interp_method
      integer :: i,j,k, orient=0, repa_mass=0, orient_mass=0, me, norient
!
      real*8, allocatable :: rin_all(:,:,:), dr(:,:,:), rout_all(:,:,:)
      real*8, allocatable :: rr(:),rr_out(:),ds(:),s(:),t(:),rrpp(:)
      real*8 :: dum
!
      character(len=80), allocatable :: fname_cor_in(:), fname_cor_out(:)
      character(len=80) :: name_cor_in, name_cor_out, dummy
      character(len=20) :: methods(4), method, form
      character(len=18) :: whoami
      character(len=8) :: keyword
!
      logical :: qprint
!
 integer :: isele, i__, iend; integer, pointer::iselct(:);character(LEN=20)::word__
      integer :: fmt, natom
      real*8, pointer :: orient_weights(:), weight(:)
      integer :: ierror
      real*8 :: u(3,3) ! rotation matrix
!
!
      interface ! to linear interpolation routine
       subroutine linear_interp(xin,yin,nin,xout,yout,nout,dydxout)
       use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
       implicit none
       integer :: nin, nout
       real*8 :: xin(nin), yin(nin), xout(nout), yout(nout)
       real*8, optional :: dydxout(nout) ! tangent computation
       real*8 :: dydx(nout)
       end subroutine linear_interp
      end interface
!
      data whoami /' SM0K_INTERPOLATE>'/
      data methods/ 'LINEAR','CUBIC SPLINE','B-SPLINE','DST'/
!
      if (.not.sm0k_initialized) call sm0k_init()
!
      qprint=MPI_COMM_STRNG.ne.MPI_COMM_NULL
      qprint=qprint.and.(ME_STRNG.eq.0)
! get interpolation specifications
! interpolation type
!
      interp_method=0
      method=get_remove_parameter(comlyn,'METH',comlen)
      length=len(method)
      length=min(max(0,length),len(method));method(length+1:)='';call adjustleft(method,(/' ',tab/));length=len_trim(method)
      if (length.ge.4) then
       if (( method(1:4).eq.'LINE'(1:4) )) then
        interp_method=linear
       elseif (( method(1:4).eq.'BSPL'(1:4) )) then
        interp_method=bspline
       elseif (( method(1:4).eq.'SPLI'(1:4) )) then
        interp_method=spline
       endif
      endif
! print summary
      if (qprint) then
       if (interp_method.gt.0) then
        length=len(methods(interp_method))
length=min(max(0,length),len(methods(interp_method)));methods(interp_method)(length+1:)='';call adjustleft(methods(interp_method),(/' ',tab/));length=len_trim(methods(interp_method))
        write(msg___,6770) whoami, methods(interp_method)(1:length) ; call plainmessage(msg___)
 6770 format(/A,' WILL INTERPOLATE USING ',A,' INTERPOLATION')
       else
        if (length.gt.0) then
         write(msg___,6771) whoami, method(1:length), whoami ; call plainmessage(msg___)
 6771 format(/A,' UNRECOGNIZED INTERPOLATION METHOD: ',A,'.',/, &
     & A, ' WILL INTERPOLATE USING LINEAR INTERPOLATION')
        else
         write(msg___,6772) whoami, whoami ; call plainmessage(msg___)
 6772 format(/A,' UNSPECIFIED INTERPOLATION METHOD.',/, &
     & A, ' WILL INTERPOLATE USING LINEAR INTERPOLATION')
        endif ! length
       endif ! interp_method
      endif ! qprint
      if (interp_method.eq.0) interp_method=linear ! choose linear interpolation as default
! process other options ccccccccccccccccccccccccccccccccccccccccccccccc
! number of input replicas
      if (find_tag(comlyn, 'NIN', comlen).gt.0) then
       num_rep_in=atoi(get_remove_parameter(comlyn, 'NIN', comlen), 0)
       if (num_rep_in.le.0) then
        if (qprint) then ; write(msg___, 6781) whoami ; call plainmessage(msg___) ; endif
 6781 format(A,' NUMBER OF INPUT REPLICAS MUST BE > 0. ',/, &
     & ' NOTHING DONE.')
        return
       else
        if (qprint) then
          write(msg___,6783) whoami, num_rep_in ; call plainmessage(msg___)
        endif
 6783 format(A,' INITIAL STRING RESOLUTION: ', I5, ' REPLICAS.')
       endif ! num_rep_in<=0
      else
        if (qprint) then ; write(msg___, 6784) whoami ; call plainmessage(msg___) ; endif
 6784 format(A,' NUMBER OF INPUT REPLICAS UNSPECIFIED',/, &
     & A,' NOTHING DONE.')
         return
      endif ! indx('NIN')
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! number of output replicas
      if (find_tag(comlyn, 'NOUT', comlen).gt.0) then
       num_rep_out=atoi(get_remove_parameter(comlyn, 'NOUT', comlen), 0)
       if (num_rep_out.le.0) then
        if (qprint) then ; write(msg___, 6782) whoami ; call plainmessage(msg___) ; endif
 6782 format(A,' NUMBER OF OUTPUT REPLICAS MUST BE > 0. ',/, &
     & ' NOTHING DONE.')
        return
       else
        if (qprint) then
          write(msg___,6785) whoami, num_rep_out ; call plainmessage(msg___)
        endif
 6785 format(A,' INTERPOLATED STRING RESOLUTION: ', I5, ' REPLICAS.')
       endif ! num_rep_in<=0
      else
        if (qprint) then ; write(msg___, 6786) whoami ; call plainmessage(msg___) ; endif
 6786 format(A,' NUMBER OF OUTPUT REPLICAS UNSPECIFIED',/, &
     & A,' NOTHING DONE.')
         return
      endif ! indx('NIN')
!ccccc
! get input coordinate file info
      name_cor_in=get_remove_parameter(comlyn,'CRIN',comlen); len_cor_in=len_trim(name_cor_in)
      if (len_cor_in.le.0) then
       if (qprint) then ; write(msg___, 6789) whoami ; call plainmessage(msg___) ; endif
 6789 format(A,' INPUT COORDINATE FILE NAME UNSPECIFIED.',/, &
     & '  NOTHING DONE.')
       return
      endif
! get output coordinate file info
      name_cor_out=get_remove_parameter(comlyn,'CROUT',comlen); len_cor_out=len_trim(name_cor_out)
      if (len_cor_out.le.0) then
       if (qprint) then ; write(msg___, 6790) whoami ; call plainmessage(msg___) ; endif
 6790 format(A,' OUTPUT COORDINATE FILE NAME UNSPECIFIED.',/, &
     & ' NOTHING DONE.')
       return
      endif ! len_cor_out
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! parse file format spec. (same for both input/output)c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (remove_tag(comlyn,'PDB',comlen).gt.0) then
        fmt=pdb
        form='FORMATTED'
      elseif ( (remove_tag(comlyn,'CARD',comlen).gt.0).or. &
     & (remove_tag(comlyn,'FORM',comlen).gt.0)) then
        fmt=charmm
        form='FORMATTED'
      else ! default
        fmt=charmm
        form='FORMATTED'
      endif
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! write summary ( same code as in SMCV interpolation )
      if (qprint) then ! note: using qprint as a root node flag, too (should change this)
!ccccc get coordinate file names
       ifile=-1 ! a valid unit number will be assigned by __OPEN_FILE
       ofile=-1
       call files_open(ifile, name_cor_in(1:len_cor_in), 'FORMATTED', 'READ')
       allocate(fname_cor_in(num_rep_in))
!
       do j=1, num_rep_in
        read(ifile,*) fname_cor_in(j)
       enddo
       call files_close(ifile)
!
       write(msg___,6791) whoami ; call plainmessage(msg___)
 6791 format(A,' COORDINATE SETS WILL BE READ FROM', &
     & ' THE FOLLOWING FILES:' )
!
       do j=1, num_rep_in
        write(msg___,'(A1,I5," ",A80)') tab, j, fname_cor_in(j) ; call plainmessage(msg___)
       enddo
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       call files_open(ofile, name_cor_out(1:len_cor_out), 'FORMATTED', 'READ')
!
       allocate(fname_cor_out(num_rep_out))
!
       do j=1, num_rep_out
        read(ifile,'(A80)') fname_cor_out(j)
       enddo
       call files_close(ofile)
!
       write(msg___,6791) whoami ; call plainmessage(msg___)
 6793 format(A,' COORDINATE SETS WILL BE WRITTEN TO', &
     & ' THE FOLLOWING FILES:' )
!
       do j=1, num_rep_out
        write(msg___,'(A1,I5," ",A80)') tab, j, fname_cor_out(j) ; call plainmessage(msg___)
       enddo
!
      endif ! qprint
!
!cccccccccccccccccccccc parse weighting/orientation options
! check for orientation options: must appear after 'ORIE'
! process selection
      i=remove_tag(comlyn,'ORIE',comlen)
      if (i.gt.0) then ! only if the ORIE directive exists
! selection text taken from corman.src
       orient=1
       j=find_tag(comlyn, 'SELE', comlen)
       if (j.gt.0.and.j.lt.i) then ! sele occurs before orie
        call warning(whoami, 'ATOM SELECTION MUST BE SPECIFIED AFTER ORIE.', 0)
         return
       endif
!
          nullify(iselct)
          isele=find_tag(comlyn, 'SELE', comlen)
          msg___=comlyn(isele:comlen) ! part of string that begins with the selection
          i__=comlen-isele+1
          iend=find_tag(msg___, 'END', i__) ! location of selection termination
          iend=iend-1+isele ! index into comlyn starting from 1
          msg___=comlyn(isele:iend) ! part of string that begins with the selection and ends before ' END'
          word__=pop_string(msg___) ! remove first word (which we know is 'SELE*') from msg___
! process selection:
          nullify(iselct)
          iselct=>system_getind(msg___)
! remove selection string from command line:
          msg___=comlyn(iend:comlen) ! command line starting with 'END' (see above)
          word__=pop_string(msg___) ! remove 'END*' e.g. 'ENDING' is ok too
          comlyn(isele:isele)=' ';
          comlyn(isele+1:)=msg___ ! selection has been removed from command
          comlen=len_trim(comlyn)
       if (associated(iselct)) then ; norient=size(iselct) ; else ; norient=0 ; endif
!
       if (norient.eq.0) then
        call warning(whoami, 'NO ATOMS SELECTED FOR ORIENTATION. WILL NOT ORIENT.', 0)
        orient=0 ! invalid or missing selection for orientation
       elseif (norient.lt.3) then
        call warning(whoami, 'FEWER THAN FOUR ATOMS SELECTED FOR ORIENTATION. WILL NOT ORIENT.', 0)
        orient=0
       endif
!
      endif ! orie specified
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! check for mass weighting
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      j=1
      do while (j.gt.0)
       mlen=comlen
       j=remove_tag(comlyn,'MASS',comlen) ! s
       if ( (orient.eq.1).and.(j.ge.i) ) then
        orient_mass=1
       else if (j.gt.0) then
        repa_mass=1
        i=i-(mlen-comlen) ! how much the string has shrunk
       endif
      enddo
!
 natom=psf_natom()
!
! root node performs the interpolation
! however, other nodes are included for compatibility with CHARMM routines
!
      if (qprint) then ! using qprint as a root flag -- not great
!
       if (repa_mass.eq.1) then ; write(msg___,6720) whoami ; call plainmessage(msg___) ; endif
 6720 format(A,' INTERPOLATION WILL USE MASS-WEIGHTING')
!
       if (orient.eq.1) then
!
        allocate(orient_weights(natom))
!
        orient_weights=0
        if (associated(iselct)) then
         orient_weights(iselct) = orient_mass * m(iselct) + (1d0-orient_mass)
         deallocate(iselct)
        endif
! print a summary
        write(keyword,'(I8)') norient
        mlen=len(keyword)
        mlen=min(max(0,mlen),len(keyword));keyword(mlen+1:)='';call adjustleft(keyword,(/' ',tab/));mlen=len_trim(keyword)
        write(msg___,6700) whoami, keyword(1:mlen) ; call plainmessage(msg___)
 6700 format(A,' STRING WILL BE ORIENTED BASED ON ',A,' ATOMS')
        if (orient_mass.eq.1) then ; write(msg___,6710) whoami ; call plainmessage(msg___) ; endif
 6710 format(A,' ORIENTATION WILL USE MASS-WEIGHTING')
       endif ! orient == 1
!ccccccccccccccccccccccccc do work cccccccccccccccccccccccccccc
       write(msg___,6974) whoami ; call plainmessage(msg___)
 6974 format(A,' READING COORDINATE FILES')
! allocate memory for ALL replicas (may not be the best solution if
! the structure is very large; however, with 0K string, you will want
! to exclude water
!
      endif ! qprint
! coordinate arrays input & output
!
      if (allocated(rin_all)) deallocate(rin_all)
      allocate(rin_all(natom,3,num_rep_in))
      rin_all=unknownf ! initialize
!
      if (allocated(rout_all)) deallocate(rout_all)
      allocate(rout_all(natom,3,num_rep_out))
      rout_all=unknownf
!
      ifile=-1 ! __OPEN_FILE will assign a valut unit #
!
!
      do j=1, num_rep_in
!
        length=len_trim(fname_cor_in(j))
        dummy=fname_cor_in(j)(1:length)
!
        if (qprint) then ; call files_open(ifile, dummy, form, 'READ') ; endif
!
        select case(fmt)
         case(charmm) ; call ch_coor_read(ifile, rcomp)
         case(pdb) ; call pdb_read(ifile, rcomp)
        end select
        rin_all(:,:,j)=transpose(rcomp)
       if (qprint) then ; call files_close(ifile) ; endif
      enddo
!
! check for undefined coordinates
      if (any(rin_all.eq.unknownf)) then
        call warning(whoami, 'WARNING: SOME INPUT COORDINATES ARE UNDEFINED AFTER READING.', 0)
      endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! orient structures if requested
!
      if (qprint) then
        if (orient.eq.1) then
! set up pairs
! pairs = RESHAPE( (/ (i, i=1,natom), &
! & (i, i=1,natom) /), (/ 2,natom /), order=(/ 2,1 /) )
! orient x based on the previous replica
         do j=2,num_rep_in
          call RMSBestFit(rin_all(:,:,j-1), rin_all(:,:,j), orient_weights, u) ! superpose r_(j-i) onto r_j
          rin_all(:,:,j)=matmul(rin_all(:,:,j),u) ! apply transpose (=inverse) of u to r_j
! old CHARMM routine:
! call rotls1( &
! & rin_all(1,1,j-1), rin_all(1,2,j-1), rin_all(1,3,j-1), &
! & rin_all(1,1,j), rin_all(1,2,j), rin_all(1,3,j), &
! & natom,pairs,natom, &
! & orient_weights, &
! & .false.,.false.)
         enddo ! j
        endif ! orient
!
        allocate(weight(natom))
!
        if (repa_mass.eq.1) then
         weight=m(1:natom)
        else
         weight=1d0
        endif
!
! do the actual interpolation -- simple, not self-consistent
! allocate memory
        if (allocated(dr)) deallocate(dr)
        allocate(dr(natom,3,num_rep_in-1))
        if (allocated(ds)) deallocate(ds)
        allocate(ds(num_rep_in-1))
        if (allocated(s)) deallocate(s)
        allocate(s(num_rep_in))
        if (allocated(t)) deallocate(t)
        allocate(t(num_rep_out))
        if (allocated(rr)) deallocate(rr)
        allocate(rr(num_rep_in))
        if (allocated(rr_out)) deallocate(rr_out)
        allocate(rr_out(num_rep_out))
        if (allocated(rrpp)) deallocate(rrpp)
        allocate(rrpp(num_rep_in))
!
! compute arclength
        dr=rin_all(:,:,2:num_rep_in)-rin_all(:,:,1:num_rep_in-1)
        s(1)=0
        do i=1,num_rep_in-1
         ds(i)=sqrt(sum(matmul( (transpose(dr(:,:,i)))**2, weight**2)))
         s(i+1)=s(i)+ds(i)
        enddo
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! normalize arclength
        do i=1,num_rep_in
         s(i)=s(i)/s(num_rep_in)
        enddo
!ccccccccccccccccccccccccc
! create uniform array
        do i=1,num_rep_out
         t(i)=1.0d0*(i-1)/(num_rep_out-1)
        enddo
!cccccccccccccc now interpolate variables cccccc
        if (interp_method.eq.spline) then
         do i=1,natom
          do j=1,3
           rr=rin_all(i,j,:)
           call spline_cubic_set(num_rep_in,s,rr,0,0,0,0,rrpp)
           do k=1,num_rep_out
            call spline_cubic_val(num_rep_in,s,rr,rrpp,t(k), &
     & rout_all(i,j,k),dum,dum)
           enddo
          enddo
         enddo
        elseif (interp_method.eq.bspline) then
         do i=1,natom
          do j=1,3
           rr=rin_all(i,j,:)
           do k=1,num_rep_out
            call spline_b_val(num_rep_in,s,rr,t(k),rout_all(i,j,k))
           enddo
          enddo
         enddo
        elseif (interp_method.eq.linear) then
         do i=1,natom
          do j=1,3
           rr=rin_all(i,j,:)
           call linear_interp(s,rr,num_rep_in,t,rr_out,num_rep_out)
           rout_all(i,j,:)=rr_out
          enddo
         enddo
        endif ! interp_method
! check for undefined coordinates
        if (any(rout_all.eq.unknownf)) then
         call warning(whoami, 'WARNING: SOME COORDINATES ARE UNDEFINED AFTER INTERPOLATION.', 0)
        endif
!
!cccccccccccc write file cccccccccccccccccccccccc
!
!
        ofile=-1 ! open_file will assign a valut unit #
        do j=1,num_rep_out
!cccccccccccc write new coordinate file
         length=len_trim(fname_cor_out(j))
         dummy=fname_cor_out(j)(1:length)
         call files_open(ofile, dummy, form, 'WRITE')
!
         rcomp=transpose(rout_all(:,:,j))
         select case(fmt)
          case(charmm) ; call ch_coor_write(ofile, rcomp, bfactor)
          case(pdb) ; call pdb_write(ofile, rcomp, occupancy, bfactor)
         end select
         call files_close(ofile)
!
        enddo ! loop over new coordinate sets
!
        deallocate(rr, rr_out, dr, s, t, ds, rrpp)
      endif ! qprint
!
      if (allocated(rin_all)) deallocate(rin_all)
      if (allocated(rout_all)) deallocate(rout_all)
      if (allocated(fname_cor_in)) deallocate(fname_cor_in )
      if (allocated(fname_cor_out)) deallocate(fname_cor_out)
!
      if (associated(orient_weights)) deallocate(orient_weights)
      if (associated(weight)) deallocate(weight)
!
      end subroutine sm0k_interpolate
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sm0k_stat_init(comlyn, comlen)
!
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use psf
      use parser
      use system, only : r, rcomp, m, bfactor, occupancy
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use system, only : system_getind
      use constants
! use psf
!
      use mpi !!**CHARMM_ONLY**!##MPI
!ccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
 character(len=80) :: msg___
!
      CHARACTER(len=*) :: COMLYN
      INTEGER :: COMLEN
!
      character(len=80) :: rform, dform, sform, fcform, cform
      integer :: klen=0, ierror, i, j
      logical :: found
      logical :: qprint, qfree
!
 integer :: isele, i__, iend; integer, pointer::iselct(:);character(LEN=20)::word__
!
      real*8 :: d, com(3)
!
      character(len=8) :: keyword
      character(len=16) :: whoami
      data whoami/' SM0K_STAT_INIT>'/
!
! interface
! function sm0k_fixed_atoms()
! integer, pointer :: sm0k_fixed_atoms(:)
! end function sm0k_fixed_atoms
! end interface
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (.not.sm0k_initialized) call sm0k_init()
!
      qprint=(ME_STRNG.eq.0).and.(MPI_COMM_STRNG.ne.MPI_COMM_NULL)
! begin
! reset iteration counter
! did the user specify it?
      stat_iteration_counter=atoi(get_remove_parameter(comlyn, 'COUN', comlen), -1)
      stat_iteration_counter=max(stat_iteration_counter,0)
      if (stat_iteration_counter.gt.0) then
       if (qprint) then ; write(msg___,639) whoami, stat_iteration_counter ; call plainmessage(msg___) ; endif
 639 format(A,' SETTING ITERATION COUNTER TO ',I7)
      endif
!
!
      if (rmsd0_funit.ge.0) then ; call files_close(rmsd0_funit) ; endif
      rmsd0_fname=''
      output_rmsd0=.false.
      if (dsdt_funit.ge.0) then ; call files_close(dsdt_funit) ; endif
      dsdt_fname=''
      output_dsdt=.false.
!
      if (c_funit.ge.0) then ; call files_close(c_funit) ; endif
      c_fname=''
      output_curvature=.false.
!
      if (s_funit.ge.0) then ; call files_close(s_funit) ; endif
      s_fname=''
      output_arclength=.false.
!
! memory allocation
      if (allocated(fixed_s)) deallocate(fixed_s) ! flags
      if (allocated(iatom_s)) deallocate(iatom_s)
      if (allocated(iatom_free_s)) deallocate(iatom_free_s)
      if (allocated(rold_s)) deallocate(rold_s)
      if (allocated(rave_s)) deallocate(rave_s)
      if (allocated(rcurrent_s)) deallocate(rcurrent_s)
      if (allocated(rcomp_s)) deallocate(rcomp_s)
      if (allocated(statWeights)) deallocate(statWeights)
      if (allocated(rcomp_o)) deallocate(rcomp_o)
      if (allocated(rold_o)) deallocate(rold_o)
      if (allocated(rave_o)) deallocate(rave_o)
      nstat=0
      qstat_orient=.false.
!
!
! is there an atom selection ?
!ccccccccccccccccc first process the RMSD-related commands
      j=find_tag(comlyn, 'SELE', comlen)
      if (j.gt.0) then
          nullify(iselct)
          isele=find_tag(comlyn, 'SELE', comlen)
          msg___=comlyn(isele:comlen) ! part of string that begins with the selection
          i__=comlen-isele+1
          iend=find_tag(msg___, 'END', i__) ! location of selection termination
          iend=iend-1+isele ! index into comlyn starting from 1
          msg___=comlyn(isele:iend) ! part of string that begins with the selection and ends before ' END'
          word__=pop_string(msg___) ! remove first word (which we know is 'SELE*') from msg___
! process selection:
          nullify(iselct)
          iselct=>system_getind(msg___)
! remove selection string from command line:
          msg___=comlyn(iend:comlen) ! command line starting with 'END' (see above)
          word__=pop_string(msg___) ! remove 'END*' e.g. 'ENDING' is ok too
          comlyn(isele:isele)=' ';
          comlyn(isele+1:)=msg___ ! selection has been removed from command
          comlen=len_trim(comlyn)
      else
       iselct=>system_getind('ALL') ! select all atoms by default
      endif
      if (associated(iselct)) then ; nstat=size(iselct) ; else ; nstat=0 ; endif
!*************************************************************************
      if (nstat.eq.0) then
        call warning(whoami, 'NO ATOMS SELECTED FOR RMSD COMPUTATION. WILL NOT COMPUTE RMSD.', 0)
        output_rmsd0=.false.
        output_dsdt=.false.
        output_rmsd_ave=.false.
      else
! determine whether structures are to be oriented before comparison
       qstat_orient=(remove_tag(comlyn,'ORIE',comlen).gt.0)
       if (qstat_orient) then
        if (repa_initialized.and.norient.gt.0) then
         if (qprint) then ; write(msg___,638) whoami ; call plainmessage(msg___) ; endif
 638 format(A,' RMSD CALCULATIONS WILL INCLUDE ORIENTATION.')
        else
      call warning(whoami, 'REPARAMETRIZATION DISABLED OR NO ORIENTATION ATOMS FOUND. WILL NOT ORIENT.', 0)
         qstat_orient=.false.
        endif ! repa_initialized
       endif ! qstat_orient
!
!!!!!!!!!!!!!! RMSD from static structure in comp
       if (remove_tag(comlyn,'RMSD',comlen).gt.0) then ! request for RMSD
        output_rmsd0=.true.
!
        if (.not.allocated(rcurrent_s)) allocate(rcurrent_s(nstat,3))
        if (.not.allocated(rcomp_s)) allocate(rcomp_s(nstat,3))
        if (qstat_orient) then
         if (.not.allocated(rcomp_o)) allocate(rcomp_o(norient,3))
        endif ! qstat_orient
!
        rmsd0_fname=get_remove_parameter(COMLYN,'RNAM',COMLEN); rmsd0_flen=len_trim(rmsd0_fname)
        if (rmsd0_flen.eq.0) then
         call warning(whoami, 'NO RMSD FILE NAME SPECIFIED. WILL WRITE TO STDOUT.', 0)
         rmsd0_funit=fout
        else
         if (remove_tag(comlyn,'RAPP',comlen).gt.0) then ! APPEND?
           rform='APPEND'
         else
           rform='WRITE'
         endif
         if (qprint) then
          rmsd0_funit=-1
          call files_open(rmsd0_funit, rmsd0_fname, 'FORMATTED', rform)
         endif
        endif
!ccccccccccc print summary
        if (qprint) then
         if (rmsd0_flen.gt.0) then
          write(msg___,660 ) whoami,rmsd0_fname(1:rmsd0_flen)
         else
          write(msg___,661 ) whoami
         endif
         call plainmessage(msg___)
        endif
 660 format(A,' WILL WRITE STRING RMSD TO FILE ',A)
 661 format(A,' WILL WRITE STRING RMSD TO STDOUT.')
!
       endif ! RMSD
!!!!!!!!!!!!!! RMSD from structure at the previous step (zts/fts)
       if (remove_tag(comlyn,'DELS',comlen).gt.0) then
        output_dsdt=.true.
!
        if (.not.allocated(rold_s)) allocate(rold_s(nstat,3)) ! for storing "old" coords
        if (.not.allocated(rcurrent_s)) allocate(rcurrent_s(nstat,3))
        if (qstat_orient) then
         if (.not.allocated(rold_o)) allocate(rold_o(norient,3))
        endif
!
        dsdt_fname=get_remove_parameter(COMLYN,'DNAM',COMLEN); dsdt_flen=len_trim(dsdt_fname)
        if (dsdt_flen.eq.0) then
         call warning(whoami, 'NO DELS FILE NAME SPECIFIED. WILL WRITE TO STDOUT.', 0)
         dsdt_funit=fout
        else
         if (remove_tag(comlyn,'DAPP',comlen).gt.0) then ! APPEND?
           dform='APPEND'
         else
           dform='WRITE'
         endif
         if (qprint) then
          dsdt_funit=-1
          call files_open(dsdt_funit, dsdt_fname, 'FORMATTED', dform)
         endif
        endif
!ccccccccccc print summary
        if (qprint) then
         if (dsdt_flen.gt.0) then
          write(msg___,650 ) whoami,dsdt_fname(1:dsdt_flen)
         else
          write(msg___,651 ) whoami
         endif
         call plainmessage(msg___)
        endif
 650 format(A,' WILL WRITE STRING RMSD(I,I+1) TO FILE ',A)
 651 format(A,' WILL WRITE STRING RMSD(I,I+1) TO STDOUT.')
!
       endif ! rmsd0
!!!!!!!!!!!!!! RMSD from average structure (zts/fts)
       if (remove_tag(comlyn,'RMSA',comlen).gt.0) then
        output_rmsd_ave=.true.
!
        if (.not.allocated(rave_s)) allocate(rave_s(nstat,3)) ! for storing average coords
        if (.not.allocated(rold_s)) allocate(rold_s(nstat,3)) ! for storing "old" coords
        rave_s=0d0
        if (.not.allocated(rcurrent_s)) allocate(rcurrent_s(nstat,3))
        if (qstat_orient) then
         if (.not.allocated(rold_o)) allocate(rold_o(norient,3))
         if (.not.allocated(rave_o)) allocate(rave_o(norient,3))
         rave_o=0d0
        endif
!
        rmsd_ave_fname=get_remove_parameter(COMLYN,'RANM',COMLEN); rmsd_ave_flen=len_trim(rmsd_ave_fname)
        if (rmsd_ave_flen.eq.0) then
         call warning(whoami, 'NO RMSA FILE NAME SPECIFIED. WILL WRITE TO STDOUT.', 0)
         rmsd_ave_funit=fout
        else
         if (remove_tag(comlyn,'RAAP',comlen).gt.0) then ! APPEND?
           rform='APPEND'
         else
           rform='WRITE'
         endif
         if (qprint) then
          rmsd_ave_funit=-1
          call files_open(rmsd_ave_funit, rmsd_ave_fname, 'FORMATTED', rform)
         endif
        endif
!ccccccccccc print summary
        if (qprint) then
         if (rmsd_ave_flen.gt.0) then
          write(msg___,6500 ) whoami,rmsd_ave_fname(1:rmsd_ave_flen)
         else
          write(msg___,6510 ) whoami
         endif
         call plainmessage(msg___)
        endif
 6500 format(A,' WILL WRITE STRING RMSD FROM AVERAGE STRUC. TO FILE ',A)
 6510 format(A,' WILL WRITE STRING RMSD FROM AVERAGE STRUC. TO STDOUT.')
       endif ! rmsd_ave
!
       if (output_rmsd_ave.or.output_rmsd0.or.output_dsdt) then
! populate iatom array
        if (.not.allocated(fixed_s)) allocate(fixed_s(nstat)) ! psf indices
        if (.not.allocated(iatom_s)) allocate(iatom_s(nstat)) ! psf indices
        if (.not.allocated(iatom_free_s)) allocate(iatom_free_s(nstat)) ! free indices for atoms with minimizer
        nstat=0
        iatom_s=0; iatom_free_s=0;
!
! currently, no support for fixed atoms -- assume their number is zero
      if (associated(iselct)) then
       iatom_s=iselct
       iatom_free_s=iatom_s*3-2 ! gives 1,4,7...
       fixed_s=.false.
       deallocate(iselct)
      endif
!
        if (output_rmsd0) then
         do i=1, nstat
           j=iatom_s(i)
           rcomp_s(i,1)=rcomp(1,j);
           rcomp_s(i,2)=rcomp(2,j);
           rcomp_s(i,3)=rcomp(3,j);
         enddo
! orientation atoms
         if (qstat_orient) then
          do i=1, norient
           j=iatom_o(i)
           rcomp_o(i,1)=rcomp(1,j);
           rcomp_o(i,2)=rcomp(2,j);
           rcomp_o(i,3)=rcomp(3,j);
          enddo
! subtract COM
          com=matmul(transpose(rcomp_o), orientWeights)
          rcomp_o(:,1)=rcomp_o(:,1)-com(1)
          rcomp_o(:,2)=rcomp_o(:,2)-com(2)
          rcomp_o(:,3)=rcomp_o(:,3)-com(3)
!
          rcomp_s(:,1)=rcomp_s(:,1)-com(1)
          rcomp_s(:,2)=rcomp_s(:,2)-com(2)
          rcomp_s(:,3)=rcomp_s(:,3)-com(3)
!
         endif ! qstat_orient
        endif
!
        if (.not.(allocated(statWeights))) allocate(statWeights(nstat))
        statWeights=1d0
!
! use mass-weighting in RMSD computation?
!
        stat_rmsd_mass=(remove_tag(comlyn,'MASS',comlen).gt.0)
        if (stat_rmsd_mass) then
         if (qprint) then ; write(msg___,640) whoami ; call plainmessage(msg___) ; endif
         do i=1,nstat
          statWeights(i)=m(iatom_s(i))*statWeights(i)
         enddo
        endif ! stat_rmsd_mass
!
        d=sum(statWeights)
        if (abs(d).gt.errtol()) then
         d=1d0/d
         statWeights=d*statWeights
        endif
!
 640 format(A, ' WILL USE MASS WEIGHTING IN RMSD CALCULATIONS.')
       endif ! output_ ...
      endif ! nstat .eq. 0
!
!!!!!!!!!!!!!! ARCLENGTH
      if (remove_tag(comlyn,'ARCL',comlen).gt.0) then
        output_arclength=.true.
        s_fname=get_remove_parameter(COMLYN,'ANAM',COMLEN); s_flen=len_trim(s_fname)
        if (s_flen.eq.0) then
         call warning(whoami, 'STRING LENGTH FILE NAME NOT SPECIFIED. WILL WRITE TO STDOUT.', 0)
         s_funit=fout
        else
         if (remove_tag(comlyn,'AAPP',comlen).gt.0) then ! APPEND?
           sform='APPEND'
         else
           sform='WRITE'
         endif
         if (qprint) then
          s_funit=-1
          call files_open(s_funit, s_fname, 'FORMATTED', sform)
         endif
        endif
!ccccccccccc print summary
        if (qprint) then
         if (s_flen.gt.0) then
          write(msg___,652) whoami,s_fname(1:s_flen)
         else
          write(msg___,653) whoami
         endif
         call plainmessage(msg___)
        endif
 652 format(A,' WILL WRITE STRING LENGTH TO FILE ',A)
 653 format(A,' WILL WRITE STRING LENGTH TO STDOUT.')
!
      endif ! ARCLENGTH
!!!!!!!!!!!!!! CURVATURE
      if (remove_tag(comlyn,'CURV',comlen).gt.0) then
        output_curvature=.true.
        c_fname=get_remove_parameter(COMLYN,'CVNM',COMLEN); c_flen=len_trim(c_fname)
        if (c_flen.eq.0) then
         call warning(whoami, 'CURVATURE FILE NAME NOT SPECIFIED. WILL WRITE TO STDOUT.', 0)
         c_funit=fout
        else
         if (remove_tag(comlyn,'CAPP',comlen).gt.0) then ! APPEND?
           cform='APPEND'
         else
           cform='WRITE'
         endif
         if (qprint) then
          c_funit=-1
          call files_open(c_funit, c_fname, 'FORMATTED', cform)
         endif
        endif
!ccccccccccc print summary
        if (qprint) then
         if (c_flen.gt.0) then
          write(msg___,6521) whoami,c_fname(1:c_flen)
         else
          write(msg___,6531) whoami
         endif
         call plainmessage(msg___)
        endif
 6521 format(A,' WILL WRITE CURVATURE TO FILE ',A)
 6531 format(A,' WILL WRITE CURVATURE TO STDOUT.')
!
      endif ! CURVATURE
! if we got this far, we are probably OK
      stat_initialized=.true.
!
      end subroutine sm0k_stat_init
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sm0k_stat(n,var)
      use bestfit
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use psf
      use parser
      use system, only : r, rcomp, m, bfactor, occupancy
!
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use constants
      use mpi
!
      implicit none
!
      integer :: n
      real*8, optional :: var(*)
! locals
      character(len=8) :: keyword
      real*8 :: u(3,3), com(3)
!
      integer :: klen, i, j, k, ifile
!
! for output
!
      integer :: me, ierror
      character(len=80) :: fmt
      integer :: fmt_len
      character(len=8) :: dummy
      logical :: found
      logical :: qprint, qroot, qmanual
      real*8 :: rmsd0, rmsd0_all(SIZE_STRNG), dsdt, dsdt_all(SIZE_STRNG)
!
      character(len=11) :: whoami
      data whoami/' SM0K_STAT>'/
!
      qroot=MPI_COMM_STRNG.ne.MPI_COMM_NULL
      qprint=qroot.and.ME_STRNG.eq.0
      qmanual=(n.eq.0)
!cccccccccccccccccccccc begin ccccccccccccccccccccccccccccc
! check if the user has made an initialization call
      if (.not.sm0k_initialized) call sm0k_init()
!
      if (.not.stat_initialized) then
       call warning(whoami, 'NO OUTPUT OPTIONS SELECTED. NOTHING DONE', 0)
       return
      endif
!
      stat_iteration_counter=stat_iteration_counter+1
! define number format string for output
!
      if (qroot) then
       write(fmt,*) SIZE_STRNG
       fmt_len=len(fmt)
       fmt_len=min(max(0,fmt_len),len(fmt));fmt(fmt_len+1:)='';call adjustleft(fmt,(/' ',tab/));fmt_len=len_trim(fmt)
      endif
!
      me=mestring
!
!
      if (output_rmsd0.or.output_rmsd_ave.or.output_dsdt) then
       do i=1, nstat
         if (qmanual.or.fixed_s(i)) then
          j=iatom_s(i)
          rcurrent_s(i,1)=r(1,j)
          rcurrent_s(i,2)=r(2,j)
          rcurrent_s(i,3)=r(3,j)
         else
          j=iatom_free_s(i);
          rcurrent_s(i,1)=var(j);j=j+1
          rcurrent_s(i,2)=var(j);j=j+1
          rcurrent_s(i,3)=var(j)
         endif
       enddo ! nstat
! orientation atoms
       if (qstat_orient) then
        do i=1,norient
          if (qmanual.or.fixed_o(i)) then ! grab coordinates from main coordinate array
           j=iatom_o(i)
           rcurrent_o(i,1)=r(1,j)
           rcurrent_o(i,2)=r(2,j)
           rcurrent_o(i,3)=r(3,j)
          else ! grab coordinates provided by minimizer
           j=iatom_free_o(i) ! x-index (y- z- indices follow)
           rcurrent_o(i,1)=var(j);j=j+1
           rcurrent_o(i,2)=var(j);j=j+1
           rcurrent_o(i,3)=var(j)
          endif
        enddo ! norient
        com=matmul(transpose(rcurrent_o), orientWeights)
        rcurrent_o(:,1)=rcurrent_o(:,1)-com(1)
        rcurrent_o(:,2)=rcurrent_o(:,2)-com(2)
        rcurrent_o(:,3)=rcurrent_o(:,3)-com(3)
        rcurrent_s(:,1)=rcurrent_s(:,1)-com(1)
        rcurrent_s(:,2)=rcurrent_s(:,2)-com(2)
        rcurrent_s(:,3)=rcurrent_s(:,3)-com(3)
       endif ! qstat_orient
      endif ! output rmsd
!
!ccccccccccccccccccccccccccccccccccccc
      if (output_rmsd0) then
        if (qroot) then
! reference structure is in the comparison set
         if (qstat_orient) then
           call RMSBestFit(rcurrent_o,rcomp_o,orientWeights,u)
! transform current structure to overlap with reference
! (if orientation is off, u=I)
           u=transpose(u)
           rcurrent_o=matmul(rcurrent_o, u) ! need to rotate both (see below)
           rcurrent_s=matmul(rcurrent_s, u)
         endif
         rmsd0=rmsd(rcurrent_s, rcomp_s, statWeights)
!
! gather !
         call mpi_gather(rmsd0,1,MPI_REAL &
     & ,rmsd0_all,1,MPI_REAL,0, &
     & MPI_COMM_STRNG, ierror)
         if (qprint) then ! root writes
           if (rmsd0_funit.eq.fout) then
      write(rmsd0_funit,'("RMSD0> ",I5," ",'//fmt(1:fmt_len)//'F11.5)') &
     & stat_iteration_counter, &
     & (rmsd0_all(i),i=1,SIZE_STRNG)
           else
            write(rmsd0_funit,'(I5," ",'//fmt(1:fmt_len)//'F11.5)')     &
     & stat_iteration_counter, &
     & (rmsd0_all(i),i=1,SIZE_STRNG)
           endif
! flush unit: close and reopen
           call files_close(rmsd0_funit)
           call files_open(rmsd0_funit, rmsd0_fname, 'FORMATTED', 'APPEND')
! done
         endif ! qprint
        endif ! qroot
      endif
! 667 format(I5,' ',100F11.5)
! 6670 format('RMSD0> ',I5,' ',100F11.5)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_dsdt) then
        if (qroot) then
         if (repa_initialized) then ! proceed only if rold is defined
          if (qstat_orient) then
           call RMSBestFit(rcurrent_o,rold_o,orientWeights,u)
           u=transpose(u)
           rcurrent_o=matmul(rcurrent_o, u)
           rcurrent_s=matmul(rcurrent_s, u)
          endif
          dsdt=rmsd(rcurrent_s, rold_s, statWeights)
! gather !
          call mpi_gather(dsdt,1,MPI_REAL &
     & ,dsdt_all,1,MPI_REAL,0, &
     & MPI_COMM_STRNG, ierror)
          if (qprint) then ! root writes
           if (dsdt_funit.eq.fout) then
      write(dsdt_funit,'("DLEN> ",I5," ",'//fmt(1:fmt_len)//'F11.5)')   &
     & stat_iteration_counter, &
     & (dsdt_all(i),i=1,SIZE_STRNG)
           else
            write(dsdt_funit,'(I5," ",'//fmt(1:fmt_len)//'F11.5)')      &
     & stat_iteration_counter, &
     & (dsdt_all(i),i=1,SIZE_STRNG)
           endif
! flush unit: close and reopen
           call files_close(dsdt_funit)
           call files_open(dsdt_funit, dsdt_fname, 'FORMATTED', 'APPEND')
! done
          endif ! qprint
         else ! repa_initialized
          call warning(whoami, 'NO REPARAMETRIZATION OPTIONS SELECTED. SKIPPING DSDT.', 0)
         endif
        endif ! qroot
! 6680 format('DLEN> ',I5,' ',100F11.5)
! 668 format(I5,' ',100F11.5)
      endif ! dsdt
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_rmsd_ave) then
        if (qroot) then
         if (repa_initialized) then ! proceed only if rave defined
          if (qstat_orient) then
           call RMSBestFit(rcurrent_o,rave_o,orientWeights,u)
           u=transpose(u)
           rcurrent_o=matmul(rcurrent_o, u)
           rcurrent_s=matmul(rcurrent_s, u)
          endif
          dsdt=rmsd(rcurrent_s, rave_s, statWeights)
! gather !
          call mpi_gather(rmsd0,1,MPI_REAL &
     & ,rmsd0_all,1,MPI_REAL,0, &
     & MPI_COMM_STRNG, ierror)
          if (qprint) then ! root writes
           if (rmsd_ave_funit.eq.fout) then
            write(rmsd_ave_funit,'("RMSD_AVE> ",I5," ",'//              &
     & fmt(1:fmt_len)//'F11.5)')                                   &
     & stat_iteration_counter, &
     & (rmsd0_all(i),i=1,SIZE_STRNG)
           else
            write(rmsd_ave_funit,'(I5," ",'//fmt(1:fmt_len)//'F11.5)')  &
     & stat_iteration_counter, &
     & (rmsd0_all(i),i=1,SIZE_STRNG)
           endif
          endif
! 6681 format('RMSD_AVE> ',I5,' ',100F11.5)
! 6682 format(I5,' ',100F11.5)
! flush unit: close and reopen
          call files_close(rmsd_ave_funit)
          call files_open(rmsd_ave_funit, rmsd_ave_fname, 'FORMATTED', 'APPEND')
        else ! repa_initialized
          call warning(whoami, 'NO REPARAMETRIZATION OPTIONS SELECTED. SKIPPING RMSD_AVE.', 0)
        endif ! repa_initialized
       endif ! qroot
! done
      endif ! RMSD from average structure
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_arclength) then
       if (qprint) then
        if (repa_initialized) then ! proceed only if arclength defined
         if (s_funit.eq.fout) then
      write(s_funit,'("ARCL> ",I5," ",'//fmt(1:fmt_len)//'F11.5)')      &
     & stat_iteration_counter, (ds(i),i=1,SIZE_STRNG-1)
         else
          write(s_funit,'(I5," ",'//fmt(1:fmt_len)//'F11.5)')           &
     & stat_iteration_counter, (ds(i),i=1,SIZE_STRNG-1)
         endif
! flush unit: close and reopen
         call files_close(s_funit)
         call files_open(s_funit, s_fname, 'FORMATTED', 'APPEND')
! done
        else
          call warning(whoami, 'NO REPARAMETRIZATION OPTIONS SELECTED. SKIPPING ARCLENGTH.', 0)
        endif
       endif ! qprint
! 669 format(I5,' ',100F11.5)
! 6690 format('ARCL> ',I5,' ',100F11.5)
      endif ! output_arclength
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_curvature) then
       if (qprint) then
        if (repa_initialized) then ! proceed only if arclength defined
         if (c_funit.eq.fout) then
      write(c_funit,'("CURV> ",I5," ",'//fmt(1:fmt_len)//'F11.5)')      &
     & stat_iteration_counter, (curv(i),i=1,SIZE_STRNG-1)
         else
          write(c_funit,'(I5," ",'//fmt(1:fmt_len)//'F11.5)')           &
     & stat_iteration_counter, (curv(i),i=1,SIZE_STRNG-1)
         endif
! flush unit: close and reopen
         call files_close(c_funit)
         call files_open(c_funit, c_fname, 'FORMATTED', 'APPEND')
! done
        else
          call warning(whoami, 'NO REPARAMETRIZATION OPTIONS SELECTED. SKIPPING CURVATURE.', 0)
        endif
       endif ! me
! 6692 format(I5,' ',100F11.5)
! 6691 format('CURV> ',I5,' ',100F11.5)
      endif ! output_curvature
!
      end subroutine sm0k_stat
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end module sm0k
!**CHARMM_ONLY**!##ENDIF
