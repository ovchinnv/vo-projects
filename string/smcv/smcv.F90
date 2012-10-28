/*#define __WRN(__WHO,__MSG) write(0,*) 'WARNING FROM: ',__WHO,': ',__MSG*/
/*#define __PRINT(__MSG) write(0,'(A)') __MSG*/
/*#define __PRINT(__MSG) call plainmessage(__MSG)*/
/*#define __PRINTL(__MSG,__LEVEL) call plainmessage(__MSG,__LEVEL)*/
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
!CHARMM Element source/stringm/smcv.src $Revision: 1.4 $
!
! string code
!
!**CHARMM_ONLY**!##IF STRINGM
!
      SUBROUTINE smcv(COMLYN,COMLEN)
!----------------------------------------------------------------------
! command parser for the string method
!----------------------------------------------------------------------
      use sm_var
      use cv_common
      use cv_frames, only: frames_init, frames_list, frames_done, &
     & frames_read_local, frames_read_global, &
     & frames_print_local, frames_print_global, &
     & frames_calc, frames, frames_initialized, &
     & frames_reset_calculate, &
     & frames_calc_align_comp
      use cv_quaternion, only: quat_reset_calculate, quat_done
!
      use sm_config
      use smcv_master, only: smcv_fill, smcv_compute_wgt, smcv_add_hist,&
     & smcv_list, smcv_compute_M, smcv_voronoi_whereami, &
     & smcv_test_grad_fd, smcv_test_parallel, smcv_test_Minv
!
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use psf
      use parser
      use system, only : r, rcomp, m, bfactor, occupancy
      use psf
      use constants
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use constants
      use mpi
      use charmmio; use pdbio; use mol_formats
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
!
      CHARACTER(len=*) :: COMLYN
      INTEGER :: COMLEN
! local variables
      integer ivver, ivv2, iorig, ileap ! for dynamics
!
      character(len=8) :: keyword
!
!
      integer :: ifile
      integer :: i,j,n, ierror
      real*8 :: k,w,gam,step,expo_memory,zval
      integer :: num_ave_samples, irep
      integer :: c1, c2, klen, strl, delta, nskip, scol=0, totcol
      integer :: ind, all, ibeg, iend
      character(len=80) :: fname
      integer :: flen
      logical :: ok
      character(len=20) :: fixbc
      character(len=6) :: whoami
      integer*4, allocatable, dimension(:) :: temp1, temp2 ! for MPI_GRAD_TYPE (alt)
      integer :: me
!
! for interpolation
!
      character(len=20) :: methods(5), method, form
      character(len=4) :: ext
      character(len=80) :: name_cv_in, name_cv_out, name_cor_in, name_cor_out
      character(len=80) :: dummy
      character(len=80), allocatable :: fname_cor_in(:), fname_cor_out(:)
      data methods/ 'LINEAR','CUBIC SPLINE','B-SPLINE','DST','LINEAR_EXACT'/
!
      integer :: int_method, length, num_rep_in, num_rep_out, &
     & len_cv_in, len_cv_out, len_cor_in, len_cor_out, ofile
      integer :: fmt
!
      logical :: interp_cv, inte_get_coor, qcomp
      real*8, allocatable :: inte_rmsd(:,:), rtemp(:,:)
      real*8 :: voro_cut, repl_x_temp
      logical :: min_rmsd, voronoi_check_map
      integer :: which(1)
      logical :: qroot, qslave, qprint
! tests
      real*8, pointer :: fd_error(:,:) ! FD gradient computation
!
      integer(kind=MPI_ADDRESS_KIND) :: lb, extent
!
      character(len=200) :: msg___(21)=(/'','','','','','','','','','','','','','','','','','','','',''/); integer :: i_
!
!
! interface to frames_align_string routine (needed since routine is not in a module and I use optional args
      interface
       subroutine frames_align_string(x,y,z,mass,min_rmsd,ind)
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
        implicit none
        real*8 :: x(:), y(:), z(:), mass(:)
        logical, optional :: min_rmsd
        integer, optional :: ind ! frame index
       end subroutine frames_align_string
!
       subroutine frame_align_rmsd(x,y,z,mass,ind)
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
        implicit none
        real*8 :: x(:), y(:), z(:), mass(:)
        integer, optional :: ind ! frame index
       end subroutine frame_align_rmsd
!
       subroutine frame_align_voro(x,y,z,mass,ind)
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
        implicit none
        real*8 :: x(:), y(:), z(:), mass(:)
        integer, optional :: ind ! frame index
       end subroutine frame_align_voro
      end interface
!
      data whoami /' SMCV>'/
!
      keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
!
      qroot=MPI_COMM_STRNG.ne.MPI_COMM_NULL
      qslave=((MPI_COMM_LOCAL.ne.MPI_COMM_NULL).and.SIZE_LOCAL.gt.1)
      qprint=qroot.and.ME_STRNG.eq.0
!
      if (SIZE_STRNG.gt.1) call MPI_BARRIER(MPI_COMM_STRNG,ierror) !!**CHARMM_ONLY**!##MPI
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (( keyword(1:4).eq.'INIT'(1:4) )) then
        call smcv_init()
        return
      endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (.not.smcv_initialized) then
        call smcv_init()
      endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (( keyword(1:4).eq.'DONE'(1:4) )) then
        call smcv_done()
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'REPA'(1:4) )) then
       if (comlen.gt.0) then ! this is an initialization call!
        call smcv_repa_init(comlyn, comlen)
       else
        call smcv_repa()
       endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'STAT'(1:4) )) then
       if (comlen.gt.0) then ! this is an initialization call!
        call smcv_stat_init(comlyn, comlen)
       else
        call smcv_stat()
       endif ! call statistics
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'DYNA'(1:4) )) then
!ccccccccccccccc PARSE OTHER DYNAMICS OPTIONS
       voronoi_hist_on=(remove_tag(comlyn,'VORO',comlen).gt.0)
       if (voronoi_hist_on) then
        voro_cut=-999d0
        if (find_tag(comlyn, 'VCUT', comlen).gt.0) then
         voro_cut=atof(get_remove_parameter(comlyn, 'VCUT', comlen), 0d0)
         if (voro_cut.le.0d0) then
          call warning(whoami, 'VCUT MUST BE POSITIVE. NOT SET.', 0)
         else
          call cv_common_voronoi_set_cutoff(voro_cut)
         endif
        endif
        voronoi_allow_cross=(remove_tag(comlyn,'VCRS',comlen).gt.0)
        if (voronoi_allow_cross) then
         voronoi_update_freq=atoi(get_remove_parameter(comlyn, 'VCRF', comlen), 0)
         if (voronoi_update_freq.le.0) then
          call warning(whoami, 'MUST SPECIFY POSITIVE VCRF. VORONOI CELL CROSSING DISABLED.', 0)
          voronoi_allow_cross=.false.
         elseif (find_tag(comlyn, 'VINI', comlen).gt.0) then ! if vini is present
          voronoi_nocross_ini=atoi(get_remove_parameter(comlyn, 'VINI', comlen), 0) ! get it
          if (voronoi_nocross_ini.le.0) then
           call warning(whoami, 'NONPOSITIVE VINI SPECIFIED. WILL SET TO ZERO.', 0)
           voronoi_nocross_ini=0
          endif ! voronoi_nocross_ini>0
         else
          voronoi_nocross_ini=0
         endif ! voronoi_nocross_ini present
        endif ! voronoi_allow_cross
!
! initialize Voronoi data
        if (.not.cv_common_voronoi_initialized) &
     & call cv_common_voronoi_init()
! standard V. calculation case -- no crossing
        compute_whereami=.false.
        if (.not.voronoi_allow_cross) then
! create standard map (unless map is present)
         if (.not.any(cv%voronoi_map.ne.-1)) then
          cv%voronoi_map=(/ (i, i=1, nstring) /)
          compute_whereami=.true. ! will be computed by dynamc routine
         endif
        endif
!
        voronoi_check_map=(remove_tag(comlyn,'CHCK',comlen).gt.0)
!
! compute whereami
!
        if (voronoi_check_map) then
         if (qprint) then
          write(msg___, 660) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
         endif
 660 FORMAT(A,' CHECKING VORONOI MAP AGAINST CURRENT COORDINATES.')
!
         compute_whereami=.false.
         call smcv_voronoi_whereami(r(1,:),r(2,:),r(3,:),m)
!
         if (any(cv%voronoi_map.ne.-1)) then
           me=cv%voronoi_map(mestring+1)
! compare me and whereami:
           if (qroot) then
            if(SIZE_STRNG.gt.1) then
             call MPI_ALLREDUCE(me.eq.cv%voronoi_whereami, ok, &
     & 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_STRNG, ierror)
            else
             ok=me.eq.cv%voronoi_whereami
            endif
           endif ! qroot
           if (qslave) then
            call mpi_bcast(ok,1,MPI_INTEGER4,0,MPI_COMM_LOCAL,ierror)
           endif
           if (.not.ok) then
            call warning(whoami, 'VORONOI MAP INCONSISTENT WITH CURRENT COORDINATES. ABORTING.', 0)
            return
           endif ! .not. ok
         else ! voronoi map invalid (or was not read); proceed anyway using current whereami
          call warning(whoami, 'VORONOI MAP CONTAINS INVALID ENTRIES.', 0)
         endif ! voronoi_map.ne.-1
!
        else
         cv%voronoi_whereami=cv%voronoi_map(mestring+1)
        endif ! voronoi_check_map
!
       endif ! voronoi_hist_on
! reset internal interation counter for smcv_master
       olditeration=0
!
       repa_on=(remove_tag(comlyn,'REPA',comlen).gt.0)
       if (repa_on) repa_freq=atoi(get_remove_parameter(comlyn, 'REPF', comlen), 0)
!
       hist_freq=atoi(get_remove_parameter(comlyn, 'HISF', comlen), 0)
!
       stat_on=(remove_tag(comlyn,'STAT',comlen).gt.0)
       if (stat_on) stat_freq=atoi(get_remove_parameter(comlyn, 'STAF', comlen), 0)
!
       evolve_cv_on=(remove_tag(comlyn,'EVOL',comlen).gt.0)
       if (evolve_cv_on) then
        evolve_freq=atoi(get_remove_parameter(comlyn, 'EVOF', comlen), 0)
        evolve_nskip=atoi(get_remove_parameter(comlyn, 'EVOS', comlen), 0)
! express in terms of history frequency
        if (hist_freq.gt.0) evolve_nskip=evolve_nskip/hist_freq
!
        evolve_step=atof(get_remove_parameter(comlyn, 'EVST', comlen), 0.0d0) ! evolution step
! ----- types of evolution
        evolve_smooth_on=(remove_tag(comlyn,'SMOO',comlen).gt.0) ! smooth trajectory
        if (evolve_smooth_on) then
         evolve_smooth_d=atoi(get_remove_parameter(comlyn, 'EVOD', comlen), 1) ! smoothing filter
        endif
!
        evolve_expo_on=(remove_tag(comlyn,'EXPO',comlen).gt.0) ! use exponential convolution
        if (evolve_expo_on) then
         evolve_expo_mem=atof(get_remove_parameter(comlyn, 'MEMO', comlen), 0.99d0)
        endif
!
        evolve_bd_on=(remove_tag(comlyn,'BD',comlen).gt.0) ! use brownian dynamics (M not used); for T-accelerated sampling
        if (evolve_bd_on) then
! evolve step specified above (will modify this section later)
         evolve_bd_T=atof(get_remove_parameter(comlyn, 'TEMP', comlen), 0d0)
        endif
!
        evolve_aver_on=(remove_tag(comlyn,'AVER',comlen).gt.0) ! z=mean(theta)
        if (evolve_aver_on) then
         num_ave_cv_samples=atoi(get_remove_parameter(comlyn, 'NAVE', comlen), 0) ! initial number of samples in the averages
! setting this large will dampen initial fluctuations
        endif
       endif
!
       restrained_on=(remove_tag(comlyn,'RSTR',comlen).gt.0)
       if (restrained_on) then
        restrained_eq_steps=atoi(get_remove_parameter(comlyn, 'REEQ', comlen), 0)
        restrained_eq0=0
! for off-path sampling
        planar_on=(remove_tag(comlyn,'PLAN',comlen).gt.0) ! restraint applied parallel to the path
       endif
!
       unrestrained_on=(remove_tag(comlyn,'URES',comlen).gt.0)
       if (unrestrained_on) then
        unrestrained_eq_steps=atoi(get_remove_parameter(comlyn, 'UREQ', comlen), 0)
        unrestrained_eq0=0
        restrained_eq0=0
       endif
!
       repl_x_on=(remove_tag(comlyn,'REX',comlen).gt.0)
       if (repl_x_on) then
        repl_x_freq=atoi(get_remove_parameter(comlyn, 'REXF', comlen), 0)
        repl_x_temp=atof(get_remove_parameter(comlyn, 'REXT', comlen), 0d0)
!
        if (repl_x_freq.le.0) then
          call warning(whoami, 'MUST SPECIFY POSITIVE REXF. REPLICA EXCHANGE IS OFF.', 0)
          repl_x_on=.false.
!
        elseif (repl_x_temp.le.0) then
          call warning(whoami, 'MUST SPECIFY POSITIVE REXT. REPLICA EXCHANGE IS OFF.', 0)
          repl_x_on=.false.
        elseif (voronoi_hist_on) then
          call warning(whoami, 'REPLICA EXCHANGE INCOMPATIBLE WITH V. TESSELATION. REX IS OFF.', 0)
          repl_x_on=.false.
        elseif (evolve_cv_on.and.mod(repl_x_freq,evolve_freq).gt.0) then
          call warning(whoami, 'REPLICA SWAP ATTEMPT FREQ. NOT A MULTIPLE OF EVOLUTION FREQ.', 0)
        elseif (repa_on.and.mod(repl_x_freq,repa_freq).gt.0) then
          call warning(whoami, 'REPLICA SWAP ATTEMPT FREQ. NOT A MULTIPLE OF REPA. FREQ.', 0)
        else ! OK
          call cv_common_rex_set_temp(repl_x_temp)
        endif
       endif ! repl_x_on
!
       if (repa_on.or.evolve_cv_on.or.repl_x_on) then ! decrease output
         string_noprint=(remove_tag(comlyn,'NOPR',comlen).gt.0)
       endif
!--------------- DONE PARSING DYNAMICS OPTIONS -----
! print summary
!cccccccccccccccccc STRING METHOD OPTIONS cccccccccccccccccccccc
       if (qprint) then
        WRITE (msg___,'(2A)') &
     & whoami, ' STRING METHOD ENABLED.' ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
        if (evolve_cv_on) then
            WRITE (msg___,'(/,2A,/,2A,I7,A)') &
     & whoami, ' STRING EVOLUTION ENABLED.', &
     & whoami, ' WILL EVOLVE AFTER EVERY ', &
     & evolve_freq,' ITERATIONS.' ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
            WRITE (msg___,'(2A,I7,A)') &
     & whoami, ' THE FIRST', evolve_nskip, &
     & ' ITERATIONS WILL NOT CONTRIBUTE TO AVERAGES.' ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
! type of evolution
            i=0;
            if (evolve_expo_on) i=i+1
            if (evolve_bd_on) i=i+1
            if (evolve_aver_on) i=i+1
            if (evolve_smooth_on) i=i+1
!
            if (i.gt.1) then
             call warning(whoami, 'MORE THAN ONE EVOLUTION SCHEME REQUESTED. WILL USE SMCV.', 0)
             evolve_expo_on=.false.
             evolve_aver_on=.false.
             evolve_smooth_on=.false.
             evolve_bd_on=.false.
            endif
!
            if (evolve_expo_on) then
               write(msg___,671) whoami, whoami, evolve_expo_mem ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 671 format(A,' CV EVOLUTION WILL BE OF THE FORM:',/, &
     & A,' Z(N+1)=A*Z(N)+(1-A)*<THETA>, A=',F9.5,'.')
            elseif (evolve_aver_on) then
               write(msg___,6710) whoami, whoami, num_ave_cv_samples ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 6710 format(A,' CV EVOLUTION WILL BE OF THE FORM:',/, &
     &A,' Z(N+1)=AVERAGE(THETA).  INITIAL NUMBER OF SAMPLES IS ',I5,'.')
            elseif (evolve_smooth_on) then
               write(msg___,672) whoami, whoami, evolve_smooth_d ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 672 format(A,' WILL EVOLVE CV BY SMOOTHING MD TRAJECTORY',/, &
     & A,' USING FILTER WIDTH D=',F8.3)
            elseif (evolve_bd_on) then
               write(msg___,6720) whoami, whoami, evolve_step, evolve_bd_T ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 6720 format(A,' WILL EVOLVE CV USING BD ADVANCEMENT',/, &
     & A,' AT T=',F8.3,' WITH STEP=',F11.5)
            else
               write(msg___,673) whoami, whoami, evolve_step ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 673 format(A,' WILL EVOLVE CV USING SMCV ADVANCEMENT ',/, &
     & A,' WITH STEP=',F11.5)
            endif ! evolve_expo
        endif ! evolve_cv_on
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (repa_on) then
          WRITE (msg___,666) whoami, repa_freq ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 666 format(A,' WILL REPARAMETRIZE STRING AFTER EVERY ',I7, &
     & ' ITERATIONS.')
        endif ! repa_on
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (hist_freq.gt.0) then ; write(msg___,667) whoami, hist_freq ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 667 format(A,' WILL SAVE CV VALUES AFTER EVERY ', I7, ' ITERATIONS.')
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (restrained_on.and..not.unrestrained_on) then
            WRITE (msg___,'(2A)') &
     & whoami, ' WILL USE RESTRAINED DYNAMICS.' ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
        endif ! restrained
        if (unrestrained_on) then
            WRITE (msg___,'(2A)') &
     & whoami, ' WILL USE UNRESTRAINED DYNAMICS.' ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
            unrestrained_eq0=0
            WRITE (msg___,669) whoami, unrestrained_eq_steps ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 669 format(A,' WILL EQUILIBRATE UNDER CV RESTRAINTS FOR ', &
     & I7, ' STEPS.')
        endif ! unrestrained
!
        if (restrained_on.or.unrestrained_on) then
            write(msg___,665) whoami, restrained_eq_steps ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
            restrained_eq0=0
 665 format(A, ' WILL ADJUST TO NEW RESTRAINTS OVER ', &
     & I7, ' STEPS.')
        endif ! restrained
!
        if (planar_on) then
            write (msg___,'(2A)') whoami, &
     & ' WILL RESTRAIN SYSTEM IN PLANE PERPENDICULAR TO PATH.' ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
        endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (stat_on) then
            write(msg___,668) whoami, stat_freq ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 668 format(A, ' WILL OUTPUT STRING STATISTICS EVERY ', &
     & I7, ' STEPS.')
        endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (repl_x_on) then
            write(msg___,691) whoami, whoami, repl_x_freq, repl_x_temp ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 691 format(A, ' WILL ATTEMPT TO EXCHANGE NEIGHBORING REPLICAS ',/ &
     & A, ' ONCE IN EVERY ',I6,' ITERATIONS AT ',F8.3, ' K.')
        endif ! repl_x_on
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (voronoi_hist_on) then
            write(msg___,670) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 670 format(A, ' WILL COMPUTE FREE ENERGY ALONG STRING ', &
     & 'USING VORONOI TESSELLATION.' )
         if (voronoi_allow_cross) then
          write(msg___,601) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
          write(msg___,602) whoami, whoami, voronoi_update_freq ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
          if (voronoi_nocross_ini.gt.0) then
           write(msg___, 603) whoami, whoami, voronoi_nocross_ini ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
          endif
 601 format(A, ' WILL ALLOW REPLICAS TO CROSS BETWEEN V. CELLS.')
 602 format(A, ' WILL UPDATE CROSSING STATISTICS ONCE IN EVERY',/, &
     & A, I6, ' ITERATIONS.')
 603 format(A, ' WILL DISALLOW CROSSING DURING THE INITIAL ',/,A,I6, &
     & ' ITERATIONS.')
         endif
         if (voro_cut.gt.0d0) then
          write(msg___,'(2A,/2A,F11.7,A)') &
     & whoami,' STRING WILL BE RESTRICTED TO STAY WITHIN ', &
     & whoami,' THE WEIGHTED DISTANCE ',voro_cut,' OF THE CELL CENTERS.' ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
         endif
         if (restrained_on.or.unrestrained_on) then
          write(msg___,'(2A,/2A)') &
     & whoami, ' STRING DYNAMICS SHOULD BE USED WITH CAUTION', &
     & whoami, ' DURING VORONOI FE COMPUTATION.' ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
         endif
         if (evolve_cv_on.or.repa_on) then
          if (.not.voronoi_allow_cross) then
           write(msg___,'(2A,/2A)') &
     & whoami, ' STRING EVOLUTION AND REPARAMETRIZATION SHOULD', &
     & whoami, ' BE USED WITH CAUTION DURING VORONOI FE COMPUTATION.' ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
          else
           write(msg___,'(2A,/2A,/2A)') &
     & whoami, ' STRING EVOLUTION AND REPARAMETRIZATION CANNOT', &
     & whoami, ' BE USED IF VORONOI CELL CROSSING IS ALLOWED.', &
     & whoami, ' VORONOI CELL CROSSING WILL BE OFF.' ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
           voronoi_allow_cross=.false.
          endif ! voronoi_allow_cross
         endif ! evolve_cv_on
        endif ! voronoi_hist_on
       endif ! root writes
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! turn on string for dynamics
       smcv_on=.true.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'ADD '(1:4) )) then ! add CV
! call MPI_COMM_RANK(MPI_COMM_WORLD, i,ierror) ! aa
! write(600+ME_GLOBAL,*) i,ME_GLOBAL, MPI_COMM_LOCAL,
! & ME_LOCAL, SIZE_LOCAL
! call MPI_BARRIER(MPI_COMM_GLOBAL, ierror)
! stop
!
       call smcv_add(comlyn, comlen) ! this routine deals with the various CV
!
! (re)compute cv index limits (for parallelization) after each addition,
! because cv%num_cv changes
!
       if (SIZE_LOCAL.gt.0) then
        if (.not. cv_common_initialized) call cv_common_init() ! make sure cv%num_cv is defined
!
        j=ceiling(1.0d0*cv%num_cv/SIZE_LOCAL) ! max. number of CV assigned to slave node
        n=ceiling(1.0d0*cv%amap%last/SIZE_LOCAL) ! max. number of amap indices assigned to slave node
!
        do i=1,SIZE_LOCAL
         cv_send_displ(i)=min((i-1)*j,cv%num_cv-1) ! cannot exceed num_cv
         cv_send_count(i)=max(0,min(j,cv%num_cv-j*(i-1))) ! how many CV I will send to CPU i
! atom map partitioning (for parallel computation of M
!
         imap_displ(i)=min((i-1)*n,cv%amap%last-1)
         imap_count(i)=max(0,min(n,cv%amap%last-n*(i-1)))
        enddo
       endif ! SIZE
! IGNORE COMMENTS BELOW
! have to do some cheating: in order to communicate, namely, I assume that the (local) number of CV is the
! same on each processor (to that some cpus might be sending "zeros" -- cv%r that is outside of the used bounds)
! this also means that j*SIZE_LOCAL has to be less than or equal to than max_cv_common.
! basically, I only use cv_send_count(1)=j from the above array
! if (SIZE_LOCAL*j.gt.max_cv_common) then
! if (qroot) then
! if (ME_STRNG.eq.0)
! & write(msg___,'(2A,I5,1A,/,2A,I5,A)' ) whoami,
! ' CV STORAGE SPACE EXCEEDED. PARALLELIZATION REQUIRES ',
! & SIZE_LOCAL*j, ' ELEMENTS', whoami, ', BUT ONLY ', max_cv_common,
! & ' ARE ALLOCATED.'
! endif
! call wrndie(-4, whoami, ' CV STORAGE SPACE EXCEEDED.')
! endif
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
! this type will change when the auxiliary atom array is resized
! I could define the type in grad_init (only once would be required)
! but I don't want to contaminate cv_common module (it would need to know sm_config)
       if (MPI_GRAD_TYPE.ne.MPI_DATATYPE_NULL) &
     & call mpi_type_free(MPI_GRAD_TYPE, ierror)
       if (MPI_GRAD_TYPE_.ne.MPI_DATATYPE_NULL) &
     & call mpi_type_free(MPI_GRAD_TYPE_, ierror)
!
       call mpi_type_vector(6*cv%amap%last,1,cv%num_cv, & ! 6 because we are taking 3 gradients and both grad arrays (see cv_common)
     & MPI_REAL,MPI_GRAD_TYPE_,ierror)
!
! indexed version
! i=6*cv%amap%last
! allocate(temp1(i), temp2(i))
! temp1=(/ (1, j=1, i) /) ! block sizes
! temp2=(/ ( (j-1)*cv%num_cv, j=1, i) /) ! offsets from zero
! call mpi_type_indexed(i, temp1, temp2,
! & MPI_REAL, MPI_GRAD_TYPE_, ierror)
! deallocate(temp1, temp2)
!
       lb=0
       extent=sizeofreal
       call mpi_type_create_resized(MPI_GRAD_TYPE_,lb,extent, &
     & MPI_GRAD_TYPE, ierror)
       call mpi_type_commit(MPI_GRAD_TYPE, ierror)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'CLEA'(1:4) )) then ! initialize
       if (qprint) then ; write(msg___,6666) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6666 format(/A,' WILL REMOVE ALL CV.')
       call frames_done()
       call quat_done()
       call cv_common_done()
       call cv_common_init()
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccc ADDITIONAL VORONOI OPTIONS (nowhere else to put them)
      elseif (( keyword(1:4).eq.'VORO'(1:4) )) then
! get voronoi command
       keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
       if (( keyword(1:4).eq.'VMAP'(1:4) )) then
        if (remove_tag(comlyn,'CALC',comlen).gt.0) then
          if (qprint) then ; write(msg___,6010) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6010 format(A,' WILL CALCULATE VORONOI MAP FROM MAIN COORDINATES.')
          call smcv_voronoi_whereami(r(1,:),r(2,:),r(3,:),m)
! put 'whereami' into the map
          if (qroot.and.SIZE_STRNG.gt.1) then
           call MPI_ALLGATHER(cv%voronoi_whereami, 1, MPI_INTEGER, &
     & cv%voronoi_map, 1, MPI_INTEGER, MPI_COMM_STRNG, ierror)
          else
           cv%voronoi_map(mestring+1)=cv%voronoi_whereami
          endif
          if (qslave) then
          call mpi_bcast(cv%voronoi_map,nstring,MPI_INTEGER,0,MPI_COMM_LOCAL,ierror)
!
          endif
! print cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        elseif (remove_tag(comlyn,'PRIN',comlen).gt.0) then
          ifile=atoi(get_remove_parameter(comlyn, 'UNIT', comlen), -1)
          FNAME=get_remove_parameter(COMLYN,'NAME',COMLEN); FLEN=len_trim(FNAME)
! note: FNAME will be UPPER CASE
          if (qroot) then
           if (flen.GT.0) then
            if (qprint) then
             call files_open(ifile, fname, 'FORMATTED', 'WRITE')
             write(msg___,6011) whoami, fname(1:flen) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
            endif
 6011 format(A,' WRITING VORONOI MAP TO FILE ',A,'.')
            call cv_common_print_voro_map(ifile)
            if (qprint) then ; call files_close(ifile) ; endif
           else
            call warning(whoami, 'FILE NAME NOT SPECIFIED. NOTHING DONE.', 0)
           endif ! flen
          endif ! qroot
! read ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        elseif (remove_tag(comlyn,'READ',comlen).gt.0) then
          ifile=atoi(get_remove_parameter(comlyn, 'UNIT', comlen), -1)
          FNAME=get_remove_parameter(COMLYN,'NAME',COMLEN); FLEN=len_trim(FNAME)
! note: FNAME will be UPPER CASE
          if (flen.GT.0) then
            if (qprint) then
            call files_open(ifile, fname, 'FORMATTED', 'WRITE')
             write(msg___,6013) whoami, fname(1:flen) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
            endif
!
 6013 format(A,' READING VORONOI MAP FROM FILE ',A,'.')
            call cv_common_read_voro_map(ifile)
            if (qprint) then ; call files_close(ifile) ; endif
           else
            call warning(whoami, 'FILE NAME NOT SPECIFIED. NOTHING DONE.', 0)
           endif ! flen
        endif ! CALC
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       elseif (( keyword(1:4).eq.'READ'(1:4) )) then
! read "restart" file that contains (1) crossing_attempt (2) crossing_accepts (3) occupancy
         ifile=atoi(get_remove_parameter(comlyn, 'UNIT', comlen), -1)
         FNAME=get_remove_parameter(COMLYN,'NAME',COMLEN); FLEN=len_trim(FNAME)
! note: FNAME will be UPPER CASE
         if (flen.GT.0) then
          if (qprint) then
           call files_open(ifile, fname, 'FORMATTED', 'WRITE')
           write(msg___,6014) whoami, fname(1:flen) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
          endif
 6014 format(A,' READING VORONOI CROSSING DATA FROM FILE ',A,'.')
          call cv_common_read_voro_data(ifile)
          if (qprint) then ; call files_close(ifile) ; endif
         else
          call warning(whoami, 'FILE NAME NOT SPECIFIED. NOTHING DONE.', 0)
         endif ! flen
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       elseif (( keyword(1:4).eq.'PRIN'(1:4) )) then
         ifile=atoi(get_remove_parameter(comlyn, 'UNIT', comlen), -1)
         FNAME=get_remove_parameter(COMLYN,'NAME',COMLEN); FLEN=len_trim(FNAME)
! note: FNAME will be UPPER CASE
         if (flen.gt.0) then
           if (qprint) then
            call files_open(ifile, fname, 'FORMATTED', 'WRITE')
            write(msg___,6015) whoami, fname(1:flen) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 6015 format(A,' WRITING VORONOI CROSSING DATA TO FILE ',A,'.')
           endif
           call cv_common_print_voro_data(ifile)
           if (qprint) then ; call files_close(ifile) ; endif
         else
            call warning(whoami, 'FILE NAME NOT SPECIFIED. NOTHING DONE.', 0)
         endif ! flen
       endif ! VMAP
!cccccccccccccccccccc FRAMES PARSER ccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'FRAM'(1:4) )) then ! frames parser
! get frames command
       keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       if (( keyword(1:4).eq.'CLEA'(1:4) )) then ! initialize
        if (qprint) then ; write(msg___,6667) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6667 format(/A,' WILL REMOVE ALL LOCAL FRAMES.')
        call frames_done()
        call frames_init()
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       elseif (( keyword(1:4).eq.'RESE'(1:4) )) then ! initialize
        if (qprint) then ; write(msg___,6650) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6650 format(/A,' WILL FORCE RECALCULATION OF FRAME AXES.')
        call frames_reset_calculate(.true.)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       elseif (( keyword(1:4).eq.'FILL'(1:4) )) then ! compute frame axes values from current coordinates
!
        qcomp=(remove_tag(comlyn,'COMP',comlen).gt.0) ! compute CV from comp set?
!
        if (frames%num_frames.lt.1) then
         call warning(whoami, 'NO FRAMES DEFINED. NOTHING DONE.', 0)
        else
! first process special option: ALIGn
! frame axes will be calculated based on the main set, but `consistently' with the comparison set;
! specifically: the frame axes are permuted such that the rotation matrix associated with transforming one frame into another is the closest
! to the best-fit-RMSD rotation matrix
         if (remove_tag(comlyn,'ALIG',comlen).gt.0) then
          if (any(rcomp(1,:).eq.unknownf)) then
           call warning(whoami, 'MAIN X SET HAS UNDEFINED VALUES. NOTHING DONE.', 0)
          elseif (any(r(1,:).eq.unknownf)) then
           call warning(whoami, 'COMPARISON X SET HAS UNDEFINED VALUES. NOTHING DONE.', 0)
          else
           if (qcomp) then
            if (qprint) then ; write(msg___,6651) whoami, whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6651 format(/A,' WILL CALCULATE FRAME AXES FROM', &
     & ' COMPARISON COORDINATES USING', &
     & /A,' BEST-FIT ALIGNMENT WITH MAIN COORDINATES.')
            do i=1, frames%num_frames
             call frames_calc_align_comp( &
     & i,rcomp(1,:),rcomp(2,:),rcomp(3,:),r(1,:),r(2,:),r(3,:),m,.true.)
            enddo
           else ! qcomp
            if (qprint) then ; write(msg___,6652) whoami, whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6652 format(/A,' WILL CALCULATE FRAME AXES FROM', &
     & ' MAIN COORDINATES USING BEST-FIT', &
     & /A,' ALIGNMENT WITH COMPARISON COORDINATES.')
            do i=1, frames%num_frames
             call frames_calc_align_comp( &
     & i,r(1,:),r(2,:),r(3,:),rcomp(1,:),rcomp(2,:),rcomp(3,:),m,.true.)
            enddo
           endif ! qcomp
           call frames_reset_calculate(.true.) ! make sure that next time frames_calc is called we recalculate axes (to be safe)
          endif ! xcomp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! process regular read
         elseif (qcomp) then
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          if (any(rcomp(1,:).eq.unknownf)) then
           call warning(whoami, 'COMPARISON X SET HAS UNDEFINED VALUES. NOTHING DONE.', 0)
          else
           if (qprint) then ; write(msg___,6656) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6656 format(/A,' WILL CALCULATE FRAME AXES FROM COMPARISON COORDINATES.')
           do i=1, frames%num_frames
            call frames_calc(i,rcomp(1,:),rcomp(2,:),rcomp(3,:),m,.true.)
           enddo
           call frames_reset_calculate(.true.) ! make sure that next time frames_calc is called we recalculate axes (to be safe)
          endif ! xcomp.eq.unknownf
         else ! qcomp false -- use main coords
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          if (any(r(1,:).eq.unknownf)) then
           call warning(whoami, 'MAIN X SET HAS UNDEFINED VALUES. NOTHING DONE.', 0)
          else
           if (qprint) then ; write(msg___,6668) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6668 format(/A,' WILL CALCULATE FRAME AXES FROM MAIN COORDINATES.')
           do i=1, frames%num_frames
            call frames_calc(i,r(1,:),r(2,:),r(3,:),m,.true.)
           enddo
           call frames_reset_calculate(.true.) ! make sure that next time frames_calc is called we recalculate axes (to be safe)
          endif ! x.eq.unknownf
         endif ! qcomp
        endif ! num_frames < 1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       elseif (( keyword(1:4).eq.'PRIN'(1:4) )) then
! can write both a local and a global file
! local is specified with 'ALL'; global is the default
        all=remove_tag(comlyn,'ALL',comlen) ! all replicas print
! prepare file
!-----------------------------------------------------------------------------
        ifile=atoi(get_remove_parameter(comlyn, 'UNIT', comlen), -1)
        FNAME=get_remove_parameter(COMLYN,'NAME',COMLEN); FLEN=len_trim(FNAME)
! note: FNAME will be UPPER CASE
!---------------------------------- OPEN FILE --------------------------------
        if (qroot) then
         if (all.gt.0.or.qprint) then
          if (flen.gt.0) then ; call files_open(ifile, fname, 'FORMATTED', 'WRITE') ; endif
          if (ifile .eq. -1) ifile=fout ! write to output stream
         endif
!---------------------------- assume file is open, write -------------------------
         if (qprint) then ; write(msg___,6669) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6669 format(/A,' WRITING LOCAL FRAME AXES.')
         if (all.eq.0) then ; call frames_print_global(ifile) ;
         else ; call frames_print_local(ifile) ; endif
         if (all.gt.0.or.qprint) then
          if (flen.gt.0) call files_close(ifile)
         endif
        endif ! qroot
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       elseif (( keyword(1:4).eq.'READ'(1:4) )) then
! can read from both a local and a global file
! local is specified with 'ALL'; default is global
        all=remove_tag(comlyn,'ALL',comlen) ! all replicas read
! prepare file
        ifile=atoi(get_remove_parameter(comlyn, 'UNIT', comlen), -1)
        FNAME=get_remove_parameter(COMLYN,'NAME',COMLEN); FLEN=len_trim(FNAME)
! note: flen will be UPPER CASE
        if (qroot) then
!cccccccccccccccccccccccccccc OPEN FILE ccccccccccccccccccccccc
         if (all.gt.0.or.qprint) then
          if (flen.gt.0) then ; call files_open(ifile, fname, 'FORMATTED', 'READ') ; endif
         endif
         if (ifile .eq. -1) then
          Ifile=5 ! read from input file
         endif
!cccccccccccccccccc assume file is open, read ccccccccccccccccccc
         if (qprint) then ; write(msg___,6670) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6670 format(A,' READING LOCAL FRAME AXES.')
         if (all.gt.0) then ; call frames_read_local(ifile) ;
         else ; call frames_read_global(ifile) ; endif
         if (all.gt.0.or.qprint) then
          if (flen.gt.0) then ; call files_close(ifile) ; endif
         endif
        endif ! qroot
! send to slaves
        if (MPI_COMM_LOCAL.ne.MPI_COMM_NULL.and. &
     & SIZE_LOCAL.gt.1) &
     & call MPI_BCAST(frames%r(:,:,:), frames%num_frames*9, &
     & MPI_REAL, 0, MPI_COMM_LOCAL, ierror)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       elseif (( keyword(1:3).eq.'ADD'(1:3) )) then
        call smcv_frame_add(comlyn, comlen)
!
! (re)compute frame index limits (for parallelization) after each addition
!
        if (SIZE_LOCAL.gt.0) then
         if (.not.frames_initialized) call frames_init() ! make sure frames%num_frames is defined
!
         j=ceiling(1.0d0*frames%num_frames/SIZE_LOCAL) ! max. number of frames assigned to slave node
         n=ceiling(1.0d0*cv%amap%last/SIZE_LOCAL) ! max. number of amap indices assigned to slave node
!
         do i=1,SIZE_LOCAL
          fr_send_displ(i)=min((i-1)*j,frames%num_frames-1) ! cannot exceed num_cv
          fr_send_count(i)=max(0,min(j,frames%num_frames-j*(i-1))) ! how many CV I will send to CPU i
! atom map partitioning (for parallel computation of M
!
          imap_displ(i)=min((i-1)*n,cv%amap%last-1)
          imap_count(i)=max(0,min(n,cv%amap%last-n*(i-1)))
         enddo
        endif
!cccc !aa
!
! write(0,*) ME_LOCAL, fr_send_displ(ME_LOCAL+1),
! & fr_send_count(ME_LOCAL+1),frames%num_frames
! stop
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       elseif (( keyword(1:4).eq.'LIST'(1:4) )) then ! list frames
        if (qprint) then ; write(msg___,6671) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6671 format(/A,' WILL LIST LOCAL FRAMES.')
       call frames_list() ! list local frames
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       elseif (( keyword(1:4).eq.'ALIG'(1:4) )) then ! iteratively invert frame vectors (v -> -v)
! to guess the best alignment along string &
! (optional) to mimimize DIST(z,theta(x))
        min_rmsd=(remove_tag(comlyn,'RMSD',comlen).gt.0) ! look for optimal RMSD alignment : i.e. minimize DIST(z,theta(x)) ?
!
        if (qprint) then ; write(msg___,6672) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6672 format(/A,' WILL ALIGN LOCAL FRAMES.')
        if (remove_tag(comlyn,'VORO',comlen).gt.0) then
          if (qprint) then ; write(msg___,6673) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6673 format(A,' WILL FIT THE VORONOI MAP.')
          call frame_align_voro(r(1,:),r(2,:),r(3,:),m)
        else
          if (min_rmsd.and.qprint) then ; write(msg___,6674) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6674 format(A,' WILL CHECK RMSD(Z,THETA[X]).')
          call frames_align_string(r(1,:),r(2,:),r(3,:),m,min_rmsd) ! subroutine moved to sm_util to resolve dependency problems
        endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       else
            write(msg___,*)'UNRECOGNIZED SUBCOMMAND: ',keyword;call warning(whoami, msg___(1), 0)
       endif ! frames parser
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'FILL'(1:4) )) then ! fill CV values from current coordinates
!
       qcomp=(remove_tag(comlyn,'COMP',comlen).gt.0)
!
       if (qcomp) then
        if (any(rcomp(1,:).eq.unknownf)) then
         call warning(whoami, 'COMPARISON X SET HAS UNDEFINED VALUES. NOTHING DONE.', 0)
        else
         if (qprint) then ; write(msg___,6657) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6657 format(/A,' WILL OBTAIN CV VALUES FROM COMPARISON COORDINATES.')
! check for column spec
         c1=atoi(get_remove_parameter(comlyn, 'COL', comlen), -1)
         if (c1.gt.0) then
          if (qprint) then ; write(msg___,6661) whoami, c1 ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
          call smcv_fill(rcomp(1,:),rcomp(2,:),rcomp(3,:),m,c1)
          call quat_reset_calculate(.true.)
          call frames_reset_calculate(.true.)
          call cv_common_unwrap_angles(c1) ! in case they are present
         else
          if (qprint) then ; write(msg___,6662) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
          call smcv_fill(rcomp(1,:),rcomp(2,:),rcomp(3,:),m)
          call quat_reset_calculate(.true.)
          call frames_reset_calculate(.true.)
          call cv_common_unwrap_angles() ! in case they are present
         endif ! c1
        endif ! x.eq.unknownf
       else ! ~qcomp -- use main coirdinates
        if (any(r(1,:).eq.unknownf)) then
         call warning(whoami, 'MAIN X SET HAS UNDEFINED VALUES. NOTHING DONE.', 0)
        else
         if (qprint) then ; write(msg___,6660) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6660 format(/A,' WILL OBTAIN CV VALUES FROM MAIN COORDINATES.')
! check for column spec
         c1=atoi(get_remove_parameter(comlyn, 'COL', comlen), -1)
         if (c1.gt.0) then
          if (qprint) then ; write(msg___,6661) whoami, c1 ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6661 format(/A,' WILL FILL COLUMN ',I3,'.')
          call smcv_fill(r(1,:),r(2,:),r(3,:),m,c1)
          call quat_reset_calculate(.true.)
          call frames_reset_calculate(.true.)
          call cv_common_unwrap_angles(c1) ! in case they are present
         else
          if (qprint) then ; write(msg___,6662) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6662 format(/A,' WILL FILL DEFAULT COLUMN.')
          call smcv_fill(r(1,:),r(2,:),r(3,:),m)
          call quat_reset_calculate(.true.)
          call frames_reset_calculate(.true.)
          call cv_common_unwrap_angles() ! in case they are present
         endif ! c1
        endif ! x.eq.unknownf
       endif ! qcomp
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'TEST'(1:4) )) then !
       if (remove_tag(comlyn,'GRAD',comlen).gt.0) then ! finite-difference gradient test
! check fd spec
        step=atof(get_remove_parameter(comlyn, 'STEP', comlen), finite_difference_d)
        if (qprint) then ; write(msg___, 7001) whoami,whoami,step,whoami,whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 7001 format(/A,' WILL TEST GRADIENTS USING FINITE DIFFERENCES', &
     & /A,' USING DX = DY = DZ = ',F15.9,'.', &
     & /A,' MAIN COORDINATE SET MUST BE DEFINED.', &
     & /A,' WILL OVERWRITE "MAIN", "ZCUR", AND "ZOLD" CV ARRAYS')
        if (any(r(1,:).eq.unknownf)) then
         call warning(whoami, 'MAIN X SET HAS UNDEFINED VALUES. NOTHING DONE.', 0)
        else
         fd_error=>smcv_test_grad_fd(r(1,:),r(2,:),r(3,:),m,step)
         if (qprint) then
          write(msg___,7002) whoami, whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 7002 format(/A,' CV#, DGRAD_X_MAX, DGRAD_Y_MAX, DGRAD_Z_MAX', &
     & /A,' ==========================================')
          do i=1, cv%num_cv
         write(msg___,'(A," ",I5," ",3(F15.9," "))')whoami,i,fd_error(i,:) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
! write(msg___,*), i, fd_error(i,:)
          enddo
         endif ! qprint
! decide whether the test was passed
         zval=abs(maxval(fd_error))
         if (zval.lt.abs(step)*1d0) then
          write(msg___,7003) whoami, zval, whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
         else
          write(msg___,7004) whoami, zval, whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
          call warning(whoami, 'FINITE DERIVATIVE TEST FAILED.', 0)
         endif ! report test result
 7003 format(/A, ' THE MAXIMUM GRADIENT ierror IS ',F15.9,', ', &
     & /A, ' WHICH IS SMALLER THAN STEP. TEST PASSED.')
 7004 format(/A, ' THE MAXIMUM GRADIENT ierror IS ',F15.9,', ', &
     & /A, ' WHICH IS NO SMALLER THAN STEP. TEST FAILED.')
         deallocate(fd_error) ! smcv_test_grad returns a pointer to an array of abs errors
        endif
       endif ! grad
!
       if (remove_tag(comlyn,'PARA',comlen).gt.0) then ! parallel communication test
        if (qprint) then ; write(msg___, 7005) whoami,whoami,whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 7005 format(/A,' WILL COMPARE PARALLEL AND SERIAL CV COMPUTATION', &
     & /A,' MAIN COORDINATE SET MUST BE DEFINED.', &
     & /A,' WILL OVERWRITE "ZCUR" CV ARRAY')
        if (any(r(1,:).eq.unknownf)) then
         call warning(whoami, 'MAIN X SET HAS UNDEFINED VALUES. NOTHING DONE.', 0)
        else
         fd_error=>smcv_test_parallel(r(1,:),r(2,:),r(3,:),m) ! use the same array as above
         if (qprint) then
          write(msg___,7006) whoami, whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 7006 format(/A,' CV#, DCV, DGRAD_X_MAX, DGRAD_Y_MAX, DGRAD_Z_MAX', &
     & /A,' ===============================================')
          do i=1, cv%num_cv
         write(msg___,'(A," ",I5," ",4(F15.9," "))')whoami,i,fd_error(i,:) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
! write(msg___,*), i, fd_error(i,:)
          enddo
         endif ! qprint
! decide whether the test was passed
         zval=abs(maxval(fd_error))
         if (zval.lt.parallel_tolerance) then
          write(msg___,7007) whoami, zval, whoami, parallel_tolerance ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
         else
          write(msg___,7008) whoami, zval, whoami, parallel_tolerance ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
          call warning(whoami, 'PARALLEL COMPUTATION TEST FAILED.', 0)
         endif ! report test result
 7007 format(/A, ' THE MAXIMUM ierror IS ',E12.5,', ', &
     & /A, ' WHICH IS SMALLER THAN ',E12.5,'. TEST PASSED.')
 7008 format(/A, ' THE MAXIMUM ierror IS ',E12.5,', ', &
     & /A, ' WHICH IS NO SMALLER THAN ',E12.5,'. TEST FAILED.')
         deallocate(fd_error) ! smcv_test_grad returns a pointer to an array of abs errors
        endif
       endif ! para
!
       if (remove_tag(comlyn,'MINV',comlen).gt.0) then ! finite-difference gradient test
        if (qprint) then ; write(msg___, 7010) whoami,whoami,whoami, whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 7010 format(/A,' WILL COMPARE M TENSOR INVERSE COMPUTATION', &
     & /A,' USING LU DECOMPOSITION AND MULTIDIAGONAL', &
     & /A,' MATRIX INVERSION.' &
     & /A,' MAIN COORDINATE SET MUST BE DEFINED.')
        if (any(r(1,:).eq.unknownf)) then
         call warning(whoami, 'MAIN X SET HAS UNDEFINED VALUES. NOTHING DONE.', 0)
        else
         zval=smcv_test_Minv(r(1,:),r(2,:),r(3,:),m)
         if (zval.lt.parallel_tolerance) then
          write(msg___,7011) whoami, zval, whoami, parallel_tolerance ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
         else
          write(msg___,7012) whoami, zval, whoami, parallel_tolerance ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
          call warning(whoami, 'M INVERSE TEST FAILED.', 0)
         endif ! report test result
!
 7011 format(/A, ' THE MAXIMUM DIFFERENCE IS ',E12.5,', ', &
     & /A, ' WHICH IS SMALLER THAN ',E12.5,'. TEST PASSED.')
 7012 format(/A, ' THE MAXIMUM DIFFERENCE IS ',E12.5,', ', &
     & /A, ' WHICH IS NO SMALLER THAN ',E12.5,'. TEST FAILED.')
        endif
       endif
! other tests will go below this line
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! specify parallel CV calculation options
      elseif (( keyword(1:4).eq.'PARA'(1:4) )) then
       do while (comlen .gt. 1)
        keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
        select case(keyword)
         case('QUAT');
          keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
          select case(keyword)
           case('YES','ON','TRUE','T','yes','on','true','t')
            keyword='ENABLED '; calc_qt_para=.true.
            if (qprint) then ; write(msg___,7009) whoami, 'QUATERNIONS', keyword ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
           case('NO','OFF','FALSE','F','no','off','false','f')
            keyword='DISABLED' ; calc_qt_para=.false.
            if (qprint) then ; write(msg___,7009) whoami, 'QUATERNIONS', keyword ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
           case default
            call warning(whoami, 'UNKNOWN OPTION SPECIFIED FOR "QUAT"', 0)
          end select
         case('FRAM');
          keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
          select case(keyword)
           case('YES','ON','TRUE','T','yes','on','true','t')
            keyword='ENABLED '; calc_fr_para=.true.
            if (qprint) then ; write(msg___,7009) whoami, 'FRAMES', keyword ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
           case('NO','OFF','FALSE','F','no','off','false','f')
            keyword='DISABLED' ; calc_fr_para=.false.
            if (qprint) then ; write(msg___,7009) whoami, 'FRAMES', keyword ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
           case default
            call warning(whoami, 'UNKNOWN OPTION SPECIFIED FOR "FRAM"', 0)
          end select
         case('COLV');
          keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
          select case(keyword)
           case('YES','ON','TRUE','T','yes','on','true','t')
            keyword='ENABLED '; calc_cv_para=.true.
            if (qprint) then ; write(msg___,7009) whoami, 'CV', keyword ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
           case('NO','OFF','FALSE','F','no','off','false','f')
            keyword='DISABLED' ; calc_cv_para=.false.
            if (qprint) then ; write(msg___,7009) whoami, 'CV', keyword ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
           case default
            call warning(whoami, 'UNKNOWN OPTION SPECIFIED FOR "COLV"', 0)
          end select
         case('MMAT');
          keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
          select case(keyword)
           case('YES','ON','TRUE','T','yes','on','true','t')
            keyword='ENABLED '; calc_Mtensor_para=.true.
            if (qprint) then ; write(msg___,7009) whoami, 'M TENSOR', keyword ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
           case('NO','OFF','FALSE','F','no','off','false','f')
            keyword='DISABLED' ; calc_Mtensor_para=.false.
            if (qprint) then ; write(msg___,7009) whoami, 'M TENSOR', keyword ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
           case default
            call warning(whoami, 'UNKNOWN OPTION SPECIFIED FOR "MMAT"', 0)
          end select
         case('VORO');
          keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
          select case(keyword)
           case('YES','ON','TRUE','T','yes','on','true','t')
            keyword='ENABLED '; calc_voronoi_para=.true.
            if (qprint) then ; write(msg___,7009) whoami, 'VORONOI NORM', keyword ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
           case('NO','OFF','FALSE','F','no','off','false','f')
            keyword='DISABLED' ; calc_voronoi_para=.false.
            if (qprint) then ; write(msg___,7009) whoami, 'VORONOI NORM', keyword ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
           case default
            call warning(whoami, 'UNKNOWN OPTION SPECIFIED FOR "VORO"', 0)
          end select
         case default
          call warning(whoami, 'UNKNOWN OPTION SPECIFIED FOR "PARA"', 0)
        end select
       enddo ! comlen
 7009 format(/A, ' PARALLEL COMPUTATION OF ',A,' ',A)
!
      elseif (( keyword(1:4).eq.'MINV'(1:4) )) then
       keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
       select case(keyword)
        case('LU','lu')
          keyword='LU'; inverse_LU=.true.
          if (qprint) then ; write(msg___,7013) whoami, keyword ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
        case('DIAG','diag')
          keyword='MULTDIAG' ; inverse_LU=.false.
          if (qprint) then ; write(msg___,7013) whoami, keyword ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
        case default
          call warning(whoami, 'UNKNOWN MATRIX INVERSION OPTION SPECIFIED.', 0)
       end select
 7013 format(/A, ' MATRIX INVERSION WILL USE ',A,' ROUTINES.')
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:3).eq.'FIX'(1:3) )) then ! tell cv_posi about fixed "virtual" replicas
       if (qprint) then ; write(msg___, 6665) whoami, whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
       fixed_bc_0=(remove_tag(comlyn,'FIRS',comlen).gt.0)
       fixed_bc_1=(remove_tag(comlyn,'LAST',comlen).gt.0)
       if (fixed_bc_0) then
         fixbc=' '
         flen=1
       else
         fixbc=' NOT '
         flen=5
       endif
       if (qprint) then ; write(msg___,6663) whoami, fixbc(1:flen) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
       if (fixed_bc_1) then
         fixbc=' '
         flen=1
       else
         fixbc=' NOT '
         flen=5
       endif
       if (qprint) then ; write(msg___,6664) whoami, fixbc(1:flen) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
       call cv_common_set_bc(fixed_bc_0, fixed_bc_1)
 6663 format(/A,' FIRST REPLICA OF STRING WILL',A,'BE FIXED.')
 6664 format(A,' LAST REPLICA OF STRING WILL',A,'BE FIXED.'/)
 6665 format(A,' WARNING: SETTING BC REQUIRES REINITIALIZATION.',/, &
     & A,' ALL CV DATA WILL BE ERASED.')
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'PRIN'(1:4) )) then
! can write both a local and a global file
! local is specified with 'ALL'; global is the default
       all=remove_tag(comlyn,'ALL',comlen) ! all replicas print
! prepare file
!-----------------------------------------------------------------------------
       ifile=atoi(get_remove_parameter(comlyn, 'UNIT', comlen), -1)
       FNAME=get_remove_parameter(COMLYN,'NAME',COMLEN); FLEN=len_trim(FNAME)
! note: FNAME will be UPPER CASE
!---------------------------------- OPEN FILE --------------------------------
       if (qroot) then
        if (all.gt.0.or.qprint) then
         if (flen.gt.0) then ; call files_open(ifile, fname, 'FORMATTED', 'WRITE') ; endif
        endif
        if (ifile .eq. -1) ifile=fout ! write to output stream
!---------------------------- assume file is open, write -------------------------
! check for column spec
        c1=atoi(get_remove_parameter(comlyn, 'COL', comlen), -1)
        if (c1.gt.0) then
         if (qprint) then ; write(msg___,6679) whoami, c1 ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6679 format(/A,' WRITING COORDINATES FROM COLUMN ',I3)
         if (all.eq.0) then ; call cv_common_print_global(ifile, c1) ;
         else ; call cv_common_print_local(ifile,c1) ; endif
        else
         if (qprint) then ; write(msg___,6689) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6689 format(/A,' WRITING COORDINATES FROM DEFAULT COLUMN.')
         if (all.eq.0) then ; call cv_common_print_global(ifile) ;
         else ; call cv_common_print_local(ifile) ; endif
        endif ! c1
        if (all.gt.0.or.qprint) then
         if (flen.gt.0) then ; call files_close(ifile) ; endif
        endif
       endif ! qroot
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'READ'(1:4) )) then
! can read from both a local and a global file
! local is specified with 'ALL'; default is global
! can also read a specific column: specify SCOL x
       all=remove_tag(comlyn,'ALL',comlen) ! all replicas read
       if (find_tag(comlyn, 'SCOL', comlen).gt.0) then
        scol=atoi(get_remove_parameter(comlyn, 'SCOL', comlen), 0)
        if (scol.ge.1) then ! need to know total number of replicas in CV file
         if (all.gt.0) then
          call warning(whoami, 'ALL AND SCOL CANNOT BOTH BE SPECIFIED.', 0)
          return
         endif
         totcol=atoi(get_remove_parameter(comlyn, 'TCOL', comlen), 0)
         if (totcol.le.0) then
          call warning(whoami, 'MUST PROVIDE TOTAL NUMBER OF COLUMNS IN CV DATA FILE.', 0)
          return
         endif
        else ! scol.ge.1
         call warning(whoami, 'SCOL MUST BE A POSITIVE INTEGER.', 0)
         return
        endif
       endif ! SCOL present
!
! prepare file
       ifile=atoi(get_remove_parameter(comlyn, 'UNIT', comlen), -1)
       FNAME=get_remove_parameter(COMLYN,'NAME',COMLEN); FLEN=len_trim(FNAME)
! note: flen will be UPPER CASE
! check for column spec
       c1=atoi(get_remove_parameter(comlyn, 'COL', comlen), -1)
       if (qroot) then
!cccccccccccccccccccccccccccc OPEN FILE ccccccccccccccccccccccc
        if (all.gt.0.or.qprint) then
         if (flen.gt.0) then ; call files_open(ifile, fname, 'FORMATTED', 'READ') ; endif
        endif
        if(ifile .eq. -1) then
         ifile=5 ! read from input file
        endif
!cccccccccccccccccc assume file is open, read ccccccccccccccccccc
        if (c1.gt.0) then ! column spec
         if (qprint) then ; write(msg___,6699) whoami, c1 ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6699 format(A,' READING COORDINATES INTO COLUMN ',I3)
         if (all.gt.0) then ; call cv_common_read_local(ifile, c1) ;
         elseif (scol.ge.1) then
          call cv_common_read_local_from_global(ifile, totcol, scol, c1)
         else; call cv_common_read_global(ifile,c1) ; endif
        else
         if (qprint) then ; write(msg___,6709) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6709 format(A,' READING COORDINATES INTO DEFAULT COLUMN.')
         if (all.gt.0) then ; call cv_common_read_local(ifile) ;
         elseif (scol.ge.1) then
          call cv_common_read_local_from_global(ifile, totcol, scol)
         else ; call cv_common_read_global(ifile) ; endif
        endif
        if (all.gt.0.or.qprint) then
         if (flen.gt.0) then ; call files_close(ifile) ; endif
        endif
       endif ! qroot
!
! broadcast to slaves
       if (c1.lt.0) c1=main ! guess what the "default column" is
       if (MPI_COMM_LOCAL.ne.MPI_COMM_NULL.and.SIZE_LOCAL.gt.1) then
! call MPI_BCAST(cv%r(1,c1), cv%num_cv, MPI_REAL,
! & 0, MPI_COMM_LOCAL, ierror)
        call mpi_bcast(cv%r(1,c1),cv%num_cv,MPI_INTEGER8,0,MPI_COMM_LOCAL,ierror)
! broadcast BC
        if (cv_common_fixed_0_bc.eq.1) then
         call mpi_bcast(cv%r_bc_0,cv%num_cv,MPI_REAL,0,MPI_COMM_LOCAL,ierror)
        endif
        if (cv_common_fixed_1_bc.eq.1) then
         call mpi_bcast(cv%r_bc_1,cv%num_cv,MPI_REAL,0,MPI_COMM_LOCAL,ierror)
        endif
       endif
! unwrap angles if present
       call cv_common_unwrap_angles(c1)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'SWAP'(1:4) )) then ! swap two columns
! read column spec
        c1=atoi(pop_string(comlyn,comlen)) ; comlen=len_trim(comlyn)
        c2=atoi(pop_string(comlyn,comlen)) ; comlen=len_trim(comlyn)
        if (qprint) then ; write(msg___,6729) whoami, c1, c2 ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6729 format(/A,' WILL SWAP COLUMNS ',I3,' AND ',I3,' ')
        call cv_common_swap(c1,c2)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'COPY'(1:4) )) then ! copy form c1 to c2
! read column spec
        c1=atoi(pop_string(comlyn,comlen)) ; comlen=len_trim(comlyn)
        c2=atoi(pop_string(comlyn,comlen)) ; comlen=len_trim(comlyn)
        if (qprint) then ; write(msg___,6739) whoami, c1, c2 ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6739 format(/A,' WILL COPY COLUMN ',I3,' TO ',I3,' ')
        call cv_common_copy(c1,c2)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'HIST'(1:4) )) then ! parse history commands
        keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
        if (( keyword(1:3).eq.'ADD'(1:3) )) then ! save current CV values to history
         if (any(r(1,:).eq.unknownf)) then
          call warning(whoami, 'MAIN X SET HAS UNDEFINED VALUES. NOTHING DONE.', 0)
         else
          if (qprint) then ; write(msg___,674) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 674 format(A,' WILL ADD CV VALUES FROM MAIN COOR. SET INTO HISTORY.')
! last argument tells routine to add the calculated cv/derivative values to the history
          call smcv_add_hist(r(1,:),r(2,:),r(3,:),m,.true.)
         endif
        elseif (( keyword(1:4).eq.'PRIN'(1:4) )) then ! print history
! can write both a local and a global file
! local is specified with 'ALL'; global is the default
         all=remove_tag(comlyn,'ALL',comlen) ! all replicas print
! prepare file
!-----------------------------------------------------------------------------
         ifile=atoi(get_remove_parameter(comlyn, 'UNIT', comlen), -1)
         FNAME=get_remove_parameter(COMLYN,'NAME',COMLEN); FLEN=len_trim(FNAME)
! note: FNAME will be UPPER CASE
! check for other spec
         nskip=atoi(get_remove_parameter(comlyn, 'SKIP', comlen), 0) ! number of entries to skip
!---------------------------------- OPEN FILE --------------------------------
         if (qroot) then
!
          if (flen.gt.0) then ; call files_open(ifile, fname, 'FORMATTED', 'WRITE') ; endif
          if (ifile .eq. -1) ifile=fout ! write to output stream
!---------------------------- assume file is open, write -------------------------
          if (qprint) then ; write(msg___,676) whoami, nskip ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 676 format(A,' WRITING CV HISTORY. SKIPPING ',I5,' ENTRIES.')
          if (all.eq.0) then;call cv_common_print_hist_global(ifile,nskip)
          else ; call cv_common_print_hist_local(ifile,nskip) ; endif
          if (flen.gt.0) then ; call files_close(ifile) ; endif
         endif ! qroot
        elseif (( keyword(1:4).eq.'SMOO'(1:4) )) then ! smooth history
! look for other spec.
         delta=atoi(get_remove_parameter(comlyn, 'DELT', comlen), 10) ! filter width
         nskip=atoi(get_remove_parameter(comlyn, 'SKIP', comlen), 0) ! number of entries to skip
!
         if (qprint) then ; write(msg___,675) whoami, delta, nskip ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 675 format(A,' SMOOTHING CV HISTORY. FILTER WIDTH =',I5,'.', &
     & / ,' SKIPPING ',I5,' ENTRIES.')
         call cv_common_smooth_hist(delta,nskip)
        elseif (( keyword(1:4).eq.'EXPO'(1:4) )) then ! CONVOLUTION W/ EXPONENTIAL
! look for other spec.
         expo_memory=atof(get_remove_parameter(comlyn, 'MEMO', comlen), 0.99d0) ! memory in the exp. conv. kernel
         nskip=atoi(get_remove_parameter(comlyn, 'SKIP', comlen), 0) ! number of entries to skip
!
         if (qprint) then ; write(msg___,701) whoami, expo_memory, nskip ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 701 format(A,' EVOLVING CV: Z(N+1)=A*Z(N)+(1-A)*<THETA>, A=',F7.3,'.',&
     & / ,' SKIPPING ',I5,' ENTRIES.')
         call cv_common_evolve_expo(expo_memory,nskip)
        endif
! done parsing 'HIST'
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'EVOL'(1:4) )) then ! evolve string using average force (SMCV)
! look for other spec: dt and
        step=atof(get_remove_parameter(comlyn, 'STEP', comlen), 0.0d0) ! evolution step
        if (qprint) then
         if (step.eq.0.0d0) then
          call warning(whoami, 'CV EVOLUTION STEP ZERO OR UNSPECIFIED.', 0)
         endif
         write(msg___,677) whoami, step
 677 format(A,' EVOLVING CV USING AVERAGE FORCE. STEP =',F7.3,'.')
        endif ! qprint
        call cv_common_evolve_smcv(step)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:3).eq.'SET'(1:3) )) then ! modify k,w,g,dt
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! first check for global options: num_ave_samples
        if (find_tag(comlyn, 'NAVE', comlen).gt.0) then
          num_ave_samples=atoi(get_remove_parameter(comlyn, 'NAVE', comlen), -1)
          if (num_ave_samples.gt.0) then
           call cv_common_set_ave_samples(num_ave_samples)
           if (qprint) then ; write(msg___,6748) whoami, num_ave_samples ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6748 format(A,' SETTING NUMBER OF SAMPLES IN THE AVERAGE SET TO ',I7)
          else
           if (qprint) then ; write(msg___,6749) whoami, num_ave_samples; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6749 format(A,' INVALID NUMBER OF SAMPLES SPECIFIED: ',I7)
          endif
        endif
! set k parallel to path (for off-path dynamics)
        if (find_tag(comlyn, 'KPAR', comlen).gt.0) then
          k=atof(get_remove_parameter(comlyn, 'KPAR', comlen), -1d0)
          if (k.ge.0d0) then
           call cv_common_set_kpara(k)
           if (qprint) then ; write(msg___,6756) whoami, k; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6756 format(A,' SETTING PARALLEL FORCE CONSTANT TO ',F11.5)
          else
           if (qprint) then ; write(msg___,6757) whoami, k; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6757 format(A,' INVALID FORCE CONSTANT SPECIFIED: ',F11.5)
          endif
        endif ! kpara
! set k perpendicular to path (for off-path dynamics)
        if (find_tag(comlyn, 'KPRP', comlen).gt.0) then
          k=atof(get_remove_parameter(comlyn, 'KPRP', comlen), -1d0)
          if (k.ge.0d0) then
           call cv_common_set_kperp(k)
           if (qprint) then ; write(msg___,6746) whoami, k; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6746 format(A,' SETTING PERPENDICULAR FORCE CONSTANT TO ',F11.5)
          else
           if (qprint) then ; write(msg___,6747) whoami, k; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6747 format(A,' INVALID FORCE CONSTANT SPECIFIED: ',F11.5)
          endif
        endif ! kperp
!
! to set k,w,g can specify atom index, or ' ALL ' to apply to all CV
! process CV selection
        ind=atoi(get_remove_parameter(comlyn, 'IND', comlen), 0)
        all=remove_tag(comlyn,'ALL',comlen)
        if (all.gt.0) then ! will loop over all cv
         ibeg=1
         iend=cv%num_cv
         if (qprint) then ; write(msg___,6750) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6750 format(A,' ALL CV INDICES SELECTED.')
        elseif (ind.gt.0.and.ind.le.cv%num_cv) then
         ibeg=ind
         iend=ind
         if (qprint) then ; write(msg___,6751) whoami, ind ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6751 format(A,' CV INDEX ',I5,' SELECTED.')
        else ! no indices specified
         call warning(whoami, ' INVALID CV INDEX SPECIFIED', 0)
         ibeg=0
         iend=-1
        endif
!
        if (iend.gt.0) then ! skip this for invalid indices
         if (find_tag(comlyn, 'FORC', comlen).gt.0) then
          k=atof(get_remove_parameter(comlyn, 'FORC', comlen), 0.0d0)
          if (qprint) then ; write(msg___,6752) whoami, k ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6752 format(A,' WILL SET K TO ',F11.5,'.')
          do i=ibeg,iend
           call cv_common_set_k(i,k)
          enddo
         endif
!
         if (find_tag(comlyn, 'GAMM', comlen).gt.0) then
          gam=atof(get_remove_parameter(comlyn, 'GAMM', comlen), 1.0d0)
          if (qprint) then ; write(msg___,6753) whoami, gam ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6753 format(A,' WILL SET GAMMA TO ',F7.3,'.')
          do i=ibeg,iend
           call cv_common_set_g(i,gam)
          enddo
         endif
!
         if (find_tag(comlyn, 'WEIG', comlen).gt.0) then ! weighting by a real
          w=atof(get_remove_parameter(comlyn, 'WEIG', comlen), -1.0d0)
          if (qprint) then ; write(msg___,6755) whoami,w ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6755 format(A,' WILL SET WEIGHT TO ',F7.3,'.')
          do i=ibeg,iend
           call cv_common_set_w(i,w)
          enddo
         endif
!
         if (find_tag(comlyn, 'ZVAL', comlen).gt.0) then ! weighting by a real
          zval=atof(get_remove_parameter(comlyn, 'ZVAL', comlen), -1.0d0)
! check replica spec
          irep=atoi(get_remove_parameter(comlyn, 'REP', comlen), -1)
          if (irep.lt.0.or.irep.ge.nstring) then
           call warning(whoami, 'REPLICA NUMBER INVALID OR UNSPECIFIED.', 0)
          else
! check column spec
           c1=atoi(get_remove_parameter(comlyn, 'COL', comlen), -1)
           if (c1.gt.0) then
            if (qprint) then
             write(msg___,6774) whoami, irep,c1, zval
             do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
            endif
 6774 format(A,' WILL SET REPLICA ',I5,' CV VALUE IN COLUMN ', &
     & I3, ' TO ',F7.3,'.')
            if (mestring.eq.irep) then ;do i=ibeg,iend
                                       call cv_common_set_r(i,zval,c1)
                                      enddo; endif
           else
            if (qprint) then ; write(msg___,6773) whoami, irep, zval ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6773 format(A,' WILL SET REPLICA ',I5,' CV VALUE IN DEFAULT COLUMN TO '&
     & ,F7.3,'.')
            if (mestring.eq.irep) then ;do i=ibeg,iend
                                       call cv_common_set_r(i,zval)
                                      enddo; endif
           endif ! colspec
          endif ! irep
         endif ! zval
!
        endif ! iend.gt.0
! done with 'SET' parsing
!cccccccccccccccccccccccccccccccccccc M matrix cccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'MMAT'(1:4) )) then
        if (remove_tag(comlyn,'CALC',comlen).gt.0) then
          if (qprint) then ; write(msg___,6854) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6854 format(A,' COMPUTING M(X) FROM ATOMIC COORDINATES.')
          call smcv_compute_M(r(1,:),r(2,:),r(3,:),m,.true.) ! compute M and M inverse
! print
        elseif (remove_tag(comlyn,'PRIN',comlen).gt.0) then
! if running in parallel, combine partial M entries
          if (qslave) then
! call MPI_ALLREDUCE(cv%M(1,1,1),cv%M(1,1,2),max_cv_common*cv%num_cv, ! will broadcast all rows, but only num_cv columns
! & MPI_REAL, MPI_SUM, MPI_COMM_LOCAL, ierror)
          else ! qslave
! cv%M(1:cv%num_cv,1:cv%num_cv,2)=cv%M(1:cv%num_cv,1:cv%num_cv,1)
          endif ! qslave
! check for inverse spec (inverse stored in 3)
          if (remove_tag(comlyn,'INV',comlen).gt.0) then
            ind=3
            keyword='INVERSE '; klen=8
            call cv_common_compute_Minv(inverse_LU)
          else
            ind=2
            keyword=' '; klen=0
          endif
!
          ifile=atoi(get_remove_parameter(comlyn, 'UNIT', comlen), -1)
          FNAME=get_remove_parameter(COMLYN,'NAME',COMLEN); FLEN=len_trim(FNAME)
! note: FNAME will be UPPER CASE
          if (qroot) then
           if (flen.GT.0) then
            call files_open(ifile, fname, 'FORMATTED', 'WRITE')
            if (qprint) then
             write(msg___,6859) whoami, keyword(1:klen), fname(1:flen)
             do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
            endif
 6859 format(A,' WRITING M MATRIX ',A,'TO FILE ',A,'.')
            call cv_common_print_M_global(ifile, ind)
            call files_close(ifile)
           else
            call warning(whoami, 'FILE NAME NOT SPECIFIED. NOTHING DONE.', 0)
           endif ! flen
          endif ! qroot
! read cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        elseif (remove_tag(comlyn,'READ',comlen).gt.0) then
! check for inverse spec (inverse stored in 3)
          if (remove_tag(comlyn,'INV',comlen).gt.0) then
            ind=3
            keyword='INVERSE '; klen=8
          else
            ind=2
            keyword=' '; klen=0
          endif
!
          ifile=atoi(get_remove_parameter(comlyn, 'UNIT', comlen), -1)
          FNAME=get_remove_parameter(COMLYN,'NAME',COMLEN); FLEN=len_trim(FNAME)
! note: FNAME will be UPPER CASE
          if (flen.gt.0) then
            if (qprint) then
            call files_open(ifile, fname, 'FORMATTED', 'WRITE')
             write(msg___,6858) whoami, keyword(1:klen), fname(1:flen)
             do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
            endif
 6858 format(A,' READING M MATRIX ',A,'FROM FILE ',A,'.')
            call cv_common_read_M_global(ifile, ind)
            if (qprint) then ; call files_close(ifile) ; endif
           else
            call warning(whoami, 'FILE NAME NOT SPECIFIED. NOTHING DONE.', 0)
           endif ! flen
! change calculation algorithm ccccccccccccccccccccccccccccccccccccc
        elseif (remove_tag(comlyn,'FAST',comlen).gt.0) then
          keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
          select case(keyword)
           case('YES','ON','TRUE','T','yes','on','true','t')
            keyword='ENABLED '; calc_Mtensor_fast=.true.
            if (qprint) then ; write(msg___,7014) whoami, keyword ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
           case('NO','OFF','FALSE','F','no','off','false','f')
            keyword='DISABLED' ; calc_Mtensor_fast=.false.
            if (qprint) then ; write(msg___,7014) whoami, keyword ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
           case default
            call warning(whoami, 'UNKNOWN OPTION SPECIFIED FOR "FAST"', 0)
          end select
 7014 format(/A,' SPARSE MATRIX ROUTINE FOR M TENSOR COMPUTATION ',A)
        endif
!cccccccccccccccccccccccccccccccccccc CV WEIGHTS cccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'WEIG'(1:4) )) then
        if (remove_tag(comlyn,'CALC',comlen).gt.0) then
          if (qprint) then ; write(msg___,6754) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6754 format(A,' COMPUTING CV WEIGHTS FROM M(X).')
          call smcv_compute_wgt(r(1,:),r(2,:),r(3,:),m)
! print
        elseif (remove_tag(comlyn,'PRIN',comlen).gt.0) then
! process output options
          ifile=atoi(get_remove_parameter(comlyn, 'UNIT', comlen), -1)
          FNAME=get_remove_parameter(COMLYN,'NAME',COMLEN); FLEN=len_trim(FNAME)
! note: FNAME will be UPPER CASE
          if (qroot) then
           if (flen.gt.0) then ; call files_open(ifile, fname, 'FORMATTED', 'WRITE') ; endif
           if (ifile .eq. -1) then
            ifile=fout ! write to output stream
            if (qprint) then ; write(msg___,6758) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6758 format(A,' WRITING CV WEIGHTS TO STDOUT.')
           else
            if (qprint) then ; write(msg___,6759) whoami, fname(1:flen) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6759 format(A,' WRITING CV WEIGHTS TO FILE ',A,'.')
           endif
           if (qprint) call cv_common_print_wgt(ifile) ! only root node writes
           if (flen.gt.0) then ; call files_close(ifile) ; endif
          endif ! qroot
! read
        elseif (remove_tag(comlyn,'READ',comlen).gt.0) then
! prepare file
          ifile=atoi(get_remove_parameter(comlyn, 'UNIT', comlen), -1)
          FNAME=get_remove_parameter(COMLYN,'NAME',COMLEN); FLEN=len_trim(FNAME)
! note: flen will be UPPER CASE
          if (qprint) then
           if (flen.gt.0) then ; call files_open(ifile, fname, 'FORMATTED', 'READ') ; endif
           if (ifile .eq. -1) then
             ifile=5 ! read from input file
             write(msg___,6760) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 6760 format(A,' READING CV WEIGHTS FROM STDIN.')
           else
             write(msg___,6761) whoami, fname(1:flen) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 6761 format(A,' READING CV WEIGHTS FROM FILE ',A,'.')
           endif
          endif ! qprint
          call cv_common_read_wgt(ifile) ! root and slave nodes enter
          if (qprint.and.flen.gt.0) then ; call files_close(ifile) ; endif
        endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'KPAR'(1:4) )) then
       if (find_tag(comlyn, 'SET', comlen).gt.0) then
         k=atof(get_remove_parameter(comlyn, 'SET', comlen), -1d0)
         if (k.ge.0d0) then
          call cv_common_set_kpara(k)
          if (qprint) then ; write(msg___,6763) whoami, k ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6763 format(A,' SETTING PARALLEL FORCE CONSTANT TO ',F11.5)
         else
          if (qprint) then ; write(msg___,6764) whoami, k ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6764 format(A,' INVALID FORCE CONSTANT SPECIFIED: ',F11.5)
         endif
       endif ! set
!
       if (remove_tag(comlyn,'CALC',comlen).gt.0) then
        if (qprint) then ; write(msg___, 6765) whoami, whoami, whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6765 format(/A, ' COMPUTING FORCE CONSTANTS FOR RESTRAINED ',/, &
     & A, ' DYNAMICS BY SCALING KPAR WITH CV WEIGHTS.',/, &
     & A, ' OVERWRITING PREVIOUSLY DEFINED FORCE CONSTANTS.')
        call cv_common_compute_k()
       endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'TANV'(1:4) )) then ! print/read/compute tangent to path
       if (remove_tag(comlyn,'CALC',comlen).gt.0) then ! calculate
        if (qprint) then ; write(msg___,6766) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6766 format(/A,' WILL COMPUTE TANGENT TO PATH.')
        call cv_common_compute_dr()
!
       elseif (remove_tag(comlyn,'READ',comlen).gt.0) then ! read from file
! prepare file
        ifile=atoi(get_remove_parameter(comlyn, 'UNIT', comlen), -1)
        FNAME=get_remove_parameter(COMLYN,'NAME',COMLEN); FLEN=len_trim(FNAME)
! note: flen will be UPPER CASE
!cccccccccccccccccccccccccccc OPEN FILE ccccccccccccccccccccccc
        if (qprint) then
         if (flen.gt.0) then ; call files_open(ifile, fname, 'FORMATTED', 'READ') ; endif
         if(ifile .eq. -1) then
          ifile=5 ! read from input file
         endif ! ifile
!cccccccccccccccccc assume file is open, read ccccccccccccccccccc
         write(msg___,6767) whoami
 6767 format(A,' READING VECTORS TANGENT TO PATH.')
        endif ! qprint
!
        call cv_common_read_dr(ifile) ! roots and slaves
!
        if (qprint.and.flen.gt.0) then ; call files_close(ifile) ; endif
!
       elseif (remove_tag(comlyn,'PRIN',comlen).gt.0) then ! print to file
! prepare file
        ifile=atoi(get_remove_parameter(comlyn, 'UNIT', comlen), -1)
        FNAME=get_remove_parameter(COMLYN,'NAME',COMLEN); FLEN=len_trim(FNAME)
! note: FNAME will be UPPER CASE
        if (qprint) then
!---------------------------------- OPEN FILE --------------------------------
         if (flen.gt.0) then ; call files_open(ifile, fname, 'FORMATTED', 'WRITE') ; endif
         if (ifile .eq. -1) ifile=fout ! write to output stream
!---------------------------- assume file is open, write -------------------------
         write(msg___,6768) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 6768 format(/A,' WRITING VECTORS TANGENT TO PATH.')
        endif ! qprint
        if (qroot) call cv_common_print_dr(ifile)
        if (qprint.and.flen.gt.0) then ; call files_close(ifile) ; endif
       endif ! TANV
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'LIST'(1:4) )) then ! list CV
       if (qprint) then ; write(msg___,6762) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6762 format(/A,' WILL LIST CV.')
       call smcv_list() ! this routine deals with the various CV
! write(0,*) 'ME_LOCAL: ',ME_LOCAL, 'SIZE_LOCAL:', SIZE_LOCAL
! write(600+ME_LOCAL, *) cv%r(1:cv%num_cv,1:main_offset)
! close(600+ME_LOCAL)
! write(700+ME_LOCAL, *) frames%r(3,3,1:frames%num_frames)
! close(700+ME_LOCAL)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! CV interpolation
      elseif (( keyword(1:4).eq.'INTE'(1:4) )) then
! note: this is done serially
! get specifications
!ccccccccccccc should the new cv file be interpolated from the CV old file?
        interp_cv=(remove_tag(comlyn,'INTERPCV',comlen).gt.0)
        if (qprint) then
         if (interp_cv) then
          write(msg___, 6801) whoami
         else
          write(msg___, 6802) whoami
         endif
         do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
        endif
!
 6801 format(A,' WILL OBTAIN NEW CV VALUES BY INTERPOLATION.')
 6802 format(A,' NEW CV VALUES WILL BE READ FROM FILE.')
! interpolation type
        if (interp_cv) then
         int_method=0
         method=get_remove_parameter(comlyn,'METH',comlen)
         length=len(method)
         length=min(max(0,length),len(method));method(length+1:)='';call adjustleft(method,(/' ',tab/));length=len_trim(method)
         if (length.ge.4) then
           if (( method(1:4).eq.'LINE'(1:4) )) then
             int_method=linear
           elseif (( method(1:4).eq.'BSPL'(1:4) )) then
             int_method=bspline
           elseif (( method(1:4).eq.'SPLI'(1:4) )) then
             int_method=spline
           elseif (( method(1:4).eq.'LIN2'(1:4) )) then
             int_method=linear_exact
           endif
         endif
! print summary
         if (qprint) then
           if (int_method.gt.0) then
             length=len(methods(int_method))
             length=min(max(0,length),len(methods(int_method)));methods(int_method)(length+1:)='';call adjustleft(methods(int_method),(/' ',tab/));length=len_trim(methods(int_method))
             write(msg___,6770) whoami, methods(int_method)(1:length)
 6770 format(/A,' WILL INTERPOLATE CV USING ',A,' INTERPOLATION')
           else
             if (length.gt.0) then
               write(msg___,6771) whoami, method(1:length), whoami
 6771 format(/A,' UNRECOGNIZED INTERPOLATION METHOD: ',A,'.',/, &
     & A, ' WILL INTERPOLATE CV USING LINEAR INTERPOLATION')
             else
              write(msg___,6772) whoami, whoami
 6772 format(/A,' UNSPECIFIED INTERPOLATION METHOD.',/, &
     & A, ' WILL INTERPOLATE CV USING LINEAR INTERPOLATION')
             endif ! length
             do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
           endif ! int_method
         endif ! prnlev
         if (int_method.eq.0) int_method=linear ! choose linear interpolation as default
        endif ! interp_cv
! process other options ccccccccccccccccccccccccccccccccccccccccccccccc
!
        if (find_tag(comlyn, 'NIN', comlen).gt.0) then
          num_rep_in=atoi(get_remove_parameter(comlyn, 'NIN', comlen), 0)
          if (num_rep_in.le.0) then
            if (qprint) then ; write(msg___, 6781) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6781 format(A,' NUMBER OF INPUT REPLICAS MUST BE > 0. ',/, &
     & ' NOTHING DONE.')
            return
          else
            name_cv_in=get_remove_parameter(comlyn,'CVIN',comlen); len_cv_in=len_trim(name_cv_in)
            if (len_cv_in.le.0) then
              if (qprint) then ; write(msg___, 6782) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6782 format(A,' INPUT CV FILE NAME UNSPECIFIED.',/, &
     & ' NOTHING DONE.')
              return
            else
              if (qprint) then ; write(msg___,6783) &
     & whoami, num_rep_in, whoami, name_cv_in(1:len_cv_in) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6783 format(A,' INITIAL STRING RESOLUTION: ', I5, ' REPLICAS.',/, &
     & A,' INPUT CV FILE IS ', A)
            endif ! len_cv_in<=0
          endif ! num_rep_in<=0
        else
          if (qprint) then ; write(msg___, 6784) whoami, whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6784 format(A,' NUMBER OF INPUT REPLICAS UNSPECIFIED',/, &
     & A,' NOTHING DONE.')
          return
        endif ! indx('NIN')
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
! process output CV specification
        if (find_tag(comlyn, 'NOUT', comlen).gt.0) then
          num_rep_out=atoi(get_remove_parameter(comlyn, 'NOUT', comlen), 0)
          if (num_rep_out.le.0) then
            if (qprint) then ; write(msg___, 6785) whoami, whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6785 format(A,' NUMBER OF OUTPUT REPLICAS MUST BE > 0. ',/,A, &
     & ' NOTHING DONE.')
            return
          else
            name_cv_out=get_remove_parameter(comlyn,'CVOUT',comlen); len_cv_out=len_trim(name_cv_out)
            if (len_cv_out.le.0) then
              if (qprint) then ; write(msg___, 6786) whoami, whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6786 format(A,' OUTPUT CV FILE NAME UNSPECIFIED.',/,A, &
     & ' NOTHING DONE.')
              return
            else
              if (qprint) then
               write(msg___,6787) whoami, num_rep_out, whoami, name_cv_out(1:len_cv_out)
               do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
              endif
 6787 format(A,' OUTPUT STRING RESOLUTION: ', I5, ' REPLICAS.',/, &
     & A,' OUPUT CV FILE IS ', A)
            endif ! len_cv_out
          endif ! num_rep_out
        else ! num_rep_out
          if (qprint) then ; write(msg___, 6788) whoami, whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6788 format(A,' NUMBER OF OUTPUT REPLICAS UNSPECIFIED',/, &
     & A,' NOTHING DONE.')
          return
        endif ! indx('NOUT')
!ccccccccccccc coordinate file specification
        inte_get_coor=(remove_tag(comlyn,'COOR',comlen).gt.0)
        if (inte_get_coor) then ! look for input and output coordinate files
          if (qprint) then ; write(msg___,6800) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6800 format(A,' WILL GENERATE REPLICA COORDINATE SETS.')
          inte_get_coor=.true.
! will also interpolate coordinates
          name_cor_in=get_remove_parameter(comlyn,'CRIN',comlen); len_cor_in=len_trim(name_cor_in) ! text file which contains a list of file names (6/20/2011)
          if (len_cor_in.le.0) then
            if (qprint) then ; write(msg___, 6789) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6789 format(A,' INPUT COORDINATES FILE NAME UNSPECIFIED.',/, &
     & '  NOTHING DONE.')
            return
          endif
!
          name_cor_out=get_remove_parameter(comlyn,'CROUT',comlen); len_cor_out=len_trim(name_cor_out)
          if (len_cor_out.le.0) then
            if (qprint) then ; write(msg___, 6790) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6790 format(A,' OUTPUT COORDINATES FILE NAME UNSPECIFIED.',/, &
     & ' NOTHING DONE.')
            return
          endif ! len_cor_out
! parse file format spec. (same for both input/output)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! write summary
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
           write(msg___,6791) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 6791 format(A,' COORDINATE SETS WILL BE READ FROM', &
     & ' THE FOLLOWING FILES:' )
!
           do j=1, num_rep_in
            write(msg___,'(A1,I5," ",A80)') char(9), j, fname_cor_in(j) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
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
           write(msg___,6791) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 6793 format(A,' COORDINATE SETS WILL BE WRITTEN TO', &
     & ' THE FOLLOWING FILES:' )
!
           do j=1, num_rep_out
            write(msg___,'(A1,I5," ",A80)') char(9), j, fname_cor_out(j) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
           enddo
!
          endif ! qprint
        endif ! 'COOR'
!
        if (.not.(interp_cv.or.inte_get_coor)) then ! nothing to do
         write(msg___,'(A," NOTHING TO DO")') whoami
         do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
         return
        endif
!ccccccccccccccccccccccccc do work
! interpolate CV first
! compute cv weights, if needed
        if (.not.cv_common_weights_initialized) then
         call warning(whoami, 'CV WEIGHTS NOT INITIALIZED. WILL COMPUTE FROM M(X)', 0)
         call smcv_compute_wgt(r(1,:),r(2,:),r(3,:),m)
        endif
!
        if (interp_cv) then
         if (qprint) then ! only the head node does this
! open CV files
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! now prepare cv units and call interpolation routine
          ifile=-1 ! a valid unit number will be assigned by __OPEN_FILE
          ofile=-1
          call files_open(ifile, name_cv_in, 'FORMATTED', 'READ')
          call files_open(ofile, name_cv_out, 'FORMATTED', 'WRITE')
!
          call cv_common_interpolate(ifile, ofile, num_rep_in, num_rep_out,&
     & int_method)
!
          call files_close(ifile)
          call files_close(ofile)
         endif ! qprint
        endif ! interp_cv
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! process coordinate files, if requested
! this code is fairly slow; things are complicated by providing compatibility with frames, which means that
! every distance computation, essentially, required an 'optimal' rearrengement of frame vectors
! this is provided by routine `frame_align_rmsd'
        if (inte_get_coor) then
! rmsd array
         if (allocated(inte_rmsd)) deallocate(inte_rmsd)
         allocate(inte_rmsd(num_rep_out,num_rep_in))
! cv array
         if (allocated(rtemp)) deallocate(rtemp)
         allocate(rtemp(max_cv_common,num_rep_out))
!
! (re-)load new cv file and store cv in rtemp
         do j=1, num_rep_out
           if (qprint) then
            call files_open(ifile, name_cv_out, 'FORMATTED', 'READ')
           endif
           call cv_common_read_local_from_global(ifile, num_rep_out, &
     & j, comp) ! this needs to be run in parallel
           rtemp(:,j)=cv%r(:,comp)
         enddo
!
         if (qprint) then ; write(msg___,6974) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6974 format(A,' READING COORDINATE FILES')
!
         do j=1, num_rep_in
! open file
            length=len_trim(fname_cor_in(j))
            dummy=fname_cor_in(j)(1:length)
!
            if (qroot) then ; call files_open(ifile, dummy, form, 'READ') ; endif
            dummy=''
            if (mestring.eq.0) then ! need a whole group to read correctly
             select case(fmt)
              case(charmm) ; call ch_coor_read(ifile, rcomp)
              case(pdb) ; call pdb_read(ifile, rcomp)
             end select
             if (qroot) then ; call files_close(ifile) ; endif
! compute "distance" to cv values
             do i=1, num_rep_out
              cv%r(:,main)=rtemp(:,i) ! CV into main set
              if (frames_initialized) &
     & call frame_align_rmsd(rcomp(1,:), rcomp(2,:), rcomp(3,:), m) ! calculate optimal frame axes in the sense of minimal rmsd
! compute the instantaneous CV realizations
              call smcv_fill(rcomp(1,:), rcomp(2,:), rcomp(3,:), m, comp) ! works in series or parallel
! compute RMSD:
              inte_rmsd(i,j)=cv_common_rmsd(comp,main)
             enddo
            endif ! mestring
! write(600,*) inte_rmsd(:,j)
         enddo ! j=1,num_rep_in
!
! reload new cv file
         do j=1, num_rep_out
           which=minloc(inte_rmsd(j,:)) ! which index corresponds to the smallest rmsd (ds)
! open the corresponding file and save under new file name
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! open file
           length=len_trim(fname_cor_in(which(1)))
           dummy=fname_cor_in(which(1))(1:length)
!
           if (mestring.eq.0) then
            if (qroot) then ; call files_open(ifile, dummy, form, 'READ') ; endif
            select case(fmt)
             case(charmm) ; call ch_coor_read(ifile, rcomp)
             case(pdb) ; call pdb_read(ifile, rcomp)
            end select
            if (qroot) then ; call files_close(ifile) ; endif
!cccccccccccc now write the same file cccccccccccccccccccccccc
! open file
            length=len_trim(fname_cor_out(j))
            dummy=fname_cor_out(j)(1:length)
!
            if (qroot) then ; call files_open(ofile, dummy, form, 'WRITE') ; endif
!
            select case(fmt)
             case(charmm) ; call ch_coor_write(ofile, rcomp, bfactor)
             case(pdb) ; call pdb_write(ofile, rcomp, occupancy, bfactor)
            end select
            if (qroot) then ; call files_close(ofile) ; endif
!
           endif ! mestring
         enddo ! loop over new coordinate sets
         if (allocated(fname_cor_in)) deallocate(fname_cor_in )
         if (allocated(fname_cor_out)) deallocate(fname_cor_out)
        endif ! inte_get_coor
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      else
            write(msg___,*)'UNRECOGNIZED SUBCOMMAND: ',keyword;call warning(whoami, msg___(1), 0)
      endif
      end subroutine smcv
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smcv_init()
      use sm_var
      use sm_config
      use cv_common, only: cv_common_initialized, cv_common_init, main, &
     & comp, zcur, instant, runave, forces2, max_cv_common
! , only:smcv_initialized, nstring, mestring,
! & cv_send_displ,cv_send_count
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use parser
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use mpi
!
      implicit none
!
   character(len=200) :: msg___(21)=(/'','','','','','','','','','','','','','','','','','','','',''/); integer :: i_
!
!
      integer :: ierror
      logical :: qroot, qslave
      integer*4 :: temp1(3), temp2(3) ! for communication
      character(len=11) :: whoami
!
      integer(kind=MPI_ADDRESS_KIND) :: lb, extent
!
      data whoami /' SMCV_INIT>'/
!
! do a basic communicator check:
      if (ME_LOCAL.eq.0.and.ME_STRNG.eq.MPI_UNDEFINED) then
        write(msg___, 111) whoami, ME_GLOBAL, whoami
 111 FORMAT(A, ' WORLD REPLICA ',I5, ' HAS ZERO GROUP ID', &
     & /,A,' BUT INVALID STRING ID (MAY BE OK).')
        do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
      elseif (ME_STRNG.ne.MPI_UNDEFINED.and. &
     & (ME_LOCAL.ne.0.or.MPI_COMM_LOCAL.eq.MPI_COMM_NULL)) then
        write(msg___, 112) whoami, ME_GLOBAL, whoami
 112 FORMAT(A, ' WORLD REPLICA ',I5, ' HAS A VALID STRING ID', &
     & /,A,' BUT A NONZERO GROUP ID. ABORTING.')
        do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
       return
      endif
!
      qroot=ME_STRNG.ne.MPI_UNDEFINED
      qslave=ME_LOCAL.ne.MPI_UNDEFINED ! (also includes roots)
!
      if (smcv_initialized) then
       if (qroot) then
        if (ME_STRNG.eq.0) then
          write(msg___,'(2A)') &
     & whoami, ' SMCV ALREADY INITIALIZED. CALL "DONE" TO CLEAN UP.'
          do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
        endif ! ME_STRNG
       endif ! qroot
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
! broadcast string size to all slave nodes
      call mpi_bcast(nstring,1,MPI_INTEGER,0,MPI_COMM_LOCAL,ierror)
      call mpi_bcast(mestring,1,MPI_INTEGER,0,MPI_COMM_LOCAL,ierror)
!
      if (qroot) then
        if (ME_STRNG.eq.0) then
          write(msg___,'(2A,I5, A)') &
     & whoami, ' FOUND ',nstring,' REPLICAS.'
          do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
        endif
      endif
!
      smcv_initialized=.true.
!
      if (.not.cv_common_initialized) call cv_common_init()
! allocate index arrays for cv index limits (parallelization)
      allocate(cv_send_displ(SIZE_LOCAL), & ! cv
     & cv_send_count(SIZE_LOCAL))
      allocate(fr_send_displ(SIZE_LOCAL), & ! frames
     & fr_send_count(SIZE_LOCAL))
      allocate(qt_send_displ(SIZE_LOCAL), & ! quat
     & qt_send_count(SIZE_LOCAL))
      allocate(imap_displ(SIZE_LOCAL), & ! imap
     & imap_count(SIZE_LOCAL))
! initialize
      cv_send_displ=0d0
      cv_send_count=0d0
      fr_send_displ=0d0
      fr_send_count=0d0
      qt_send_displ=0d0
      qt_send_count=0d0
      imap_displ=0d0
      imap_count=0d0
!
! define/update derived MPI types:
! these special types are for communicating CV;
! I chose to index cv%r in a way that is not really
! suitable for communication/parallelization, hence the special
! types with custom blocks/strides/extents
!ccccccccccccccccc two cv values (main+comp) ccccccccccccccccc
!
! call mpi_type_indexed(2,(/1,1/),
! & (/max_cv_common*(main-1),
! & max_cv_common*(comp-1)/),
! & MPI_REAL, MPI_CV_TYPE2_, ierror)
      temp1=(/1,1,0/)
      temp2=(/max_cv_common*(main-1), &
     & max_cv_common*(comp-1),0/)
      call mpi_type_indexed(2, temp1, temp2, &
     & MPI_REAL, MPI_CV_TYPE2_, ierror)
! corresponding resized type (modified extent)
      lb=0
      extent=sizeofreal
      call mpi_type_create_resized(MPI_CV_TYPE2_,lb,extent, &
     & MPI_CV_TYPE2, ierror)
      call mpi_type_commit(MPI_CV_TYPE2, ierror)
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccc three cv values cccccccccccccc
! call mpi_type_indexed(3,(/1,1,1/),
! & (/max_cv_common*(zcur-1),
! & max_cv_common*(instant-1),
! & max_cv_common*(forces2-1)/), ! strides: i.e. take 1 element from zcur, 1 from inst., 1 from forces2
! & MPI_REAL, MPI_CV_TYPE3_, ierror)
      temp1=(/1,1,1/)
      temp2=(/max_cv_common*(zcur-1), &
     & max_cv_common*(instant-1), &
     & max_cv_common*(forces2-1)/)
      call mpi_type_indexed(3, temp1, temp2, &
     & MPI_REAL, MPI_CV_TYPE3_, ierror)
! write(0,*) temp2 !aa
! stop
! corresponding resized type
      lb=0
      extent=sizeofreal
      call mpi_type_create_resized(MPI_CV_TYPE3_,lb,extent, &
     & MPI_CV_TYPE3, ierror)
      call mpi_type_commit(MPI_CV_TYPE3, ierror)
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccc three cv values (different from above) cccccccccccccc
! call mpi_type_indexed(3,(/1,1,1/),
! & (/max_cv_common*(runave-1),
! & max_cv_common*(instant-1),
! & max_cv_common*(forces2-1)/),
! & MPI_REAL, MPI_CV_TYPE3I_, ierror)
      temp1=(/1,1,1/)
      temp2=(/max_cv_common*(runave-1), &
     & max_cv_common*(instant-1), &
     & max_cv_common*(forces2-1)/)
      call mpi_type_indexed(3, temp1, temp2, &
     & MPI_REAL, MPI_CV_TYPE3I_, ierror)
! corresponding resized type (note change of extent)
      lb=0
      extent=sizeofreal
      call mpi_type_create_resized(MPI_CV_TYPE3I_,lb,extent, &
     & MPI_CV_TYPE3I, ierror)
      call mpi_type_commit(MPI_CV_TYPE3I, ierror)
!
      MPI_GRAD_TYPE =MPI_DATATYPE_NULL
      MPI_GRAD_TYPE_=MPI_DATATYPE_NULL
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine smcv_init
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smcv_done()
      use cv_common, only: cv_common_done
      use sm_var,only: smcv_initialized, nstring, mestring
      use sm_config,only: cv_send_displ, cv_send_count, &
     & fr_send_displ,fr_send_count, &
     & imap_displ,imap_count, &
     & qt_send_displ,qt_send_count, &
     & MPI_CV_TYPE2, MPI_CV_TYPE2_, &
     & MPI_CV_TYPE3, MPI_CV_TYPE3_, &
     & MPI_CV_TYPE3I, MPI_CV_TYPE3I_, &
     & MPI_GRAD_TYPE, MPI_GRAD_TYPE_
!
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use parser
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use mpi
!
      implicit none
!
 character(len=200) :: msg___(21)=(/'','','','','','','','','','','','','','','','','','','','',''/); integer :: i_
!
      character(len=11) :: whoami
      integer :: ierror
!
      data whoami /' SMCV_DONE>'/
!
      if (MPI_COMM_STRNG.ne.MPI_COMM_NULL.and.ME_STRNG.eq.0) then
        write(msg___,'(2A,I5, A)') whoami, ' CLEANING UP.' ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
      endif
!
      call cv_common_done()
      nstring=-1
      mestring=-1
!
!
! deallocate index arrays for cv index limits (parallelization)
      if (smcv_initialized) then
       deallocate(cv_send_displ,cv_send_count, &
     & fr_send_displ,fr_send_count, &
     & imap_displ,imap_count, &
     & qt_send_displ,qt_send_count)
! free MPI_TYPES
       if (MPI_CV_TYPE2.ne.MPI_DATATYPE_NULL) &
     & call mpi_type_free(MPI_CV_TYPE2,ierror);
       if (MPI_CV_TYPE2_.ne.MPI_DATATYPE_NULL) &
     & call mpi_type_free(MPI_CV_TYPE2_,ierror);
!
       if (MPI_CV_TYPE3.ne.MPI_DATATYPE_NULL) &
     & call mpi_type_free(MPI_CV_TYPE3,ierror);
       if (MPI_CV_TYPE3_.ne.MPI_DATATYPE_NULL) &
     & call mpi_type_free(MPI_CV_TYPE3_,ierror);
!
       if (MPI_CV_TYPE3I.ne.MPI_DATATYPE_NULL) &
     & call mpi_type_free(MPI_CV_TYPE3I,ierror);
       if (MPI_CV_TYPE3I_.ne.MPI_DATATYPE_NULL) &
     & call mpi_type_free(MPI_CV_TYPE3I_,ierror);
!
       if (MPI_GRAD_TYPE.ne.MPI_DATATYPE_NULL) &
     & call mpi_type_free(MPI_GRAD_TYPE,ierror);
       if (MPI_GRAD_TYPE_.ne.MPI_DATATYPE_NULL) &
     & call mpi_type_free(MPI_GRAD_TYPE_,ierror);
      endif
!
! what else ?
!
      smcv_initialized=.false.
!
      end subroutine smcv_done
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smcv_repa_init(COMLYN, COMLEN)
! initialize string reparametrization
!
      use sm_var
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use parser
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use mpi
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
!
      CHARACTER(LEN=*) :: COMLYN
      INTEGER :: COMLEN
!
      character(len=20) :: methods(5)
      character(len=16) :: whoami
      character(len=8) :: keyword
      data methods/ 'LINEAR','CUBIC SPLINE','B-SPLINE','DST',&
     & 'LINEAR_EXACT'/
! selection array
      integer :: qlinear, qspline, qbspline, qdst, qlinear_exact
      integer :: mlen
!
! declare functions here
!
      logical :: qroot, qprint
!
 character(len=200) :: msg___(21)=(/'','','','','','','','','','','','','','','','','','','','',''/); integer :: i_
!
      data whoami /' SMCV_REPA_INIT>'/
!
      qroot=MPI_COMM_STRNG.ne.MPI_COMM_NULL
      qprint=qroot.and.ME_STRNG.eq.0
! begin
! reset variables
      qspline=0
      qbspline=0
      qlinear=0
      qdst=0
      qlinear_exact=0
      dst_cutoff=0d0
      interp_method=0
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
       if (dst_cutoff.lt.0.0d0) then
        if (qprint) then ; write(msg___,664) whoami, whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 664 FORMAT(A,' DST REQUESTED BUT FILTER CUTOFF', &
     & A, ' NOT SPECIFIED.',/,' WILL USE 0.500')
        dst_cutoff=0.5d0
       endif
      endif
      if (remove_tag(comlyn,'LIN2',comlen).gt.0) then
       qlinear_exact=1
       interp_method=linear_exact
      endif
!ccccccc CHECK FOR MULTIPLE OPTIONS
      if ((qspline+qlinear+qbspline+qdst+qlinear_exact) .eq. 0) then
       if (qprint) then ; write(msg___,665) whoami, whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 665 FORMAT(A,' INTERPOLATION METHOD NOT SPECIFIED.',/, &
     & A,' WILL USE LINEAR INTERPOLATION.')
       interp_method=linear
      elseif ((qspline+qlinear+qbspline+qdst+qlinear_exact) .gt. 1) then
       call warning(whoami, 'TOO MANY INTERPOLATION OPTIONS.', 0)
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! did the user specify a tolerance?
      if (interp_method.ne.linear_exact) then ! options below are invalid for exact interpolation
       def=atof(get_remove_parameter(comlyn, 'DEFI', comlen), 1.1d0)
       if (def.lt.1.0d0) then
         call warning(whoami, 'INTERPOLATION TOLERANCE MUST BE >= 1.', 0)
! return
       endif
! did the user specify a maximum number of iterations?
       iterations=atoi(get_remove_parameter(comlyn, 'ITER', comlen), 10)
      else
       def=0d0
       iterations=0
      endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! print summary
      if (qprint) then
       mlen=len(methods(interp_method))
       if (interp_method.eq.linear_exact) then
        write(msg___,666) whoami,methods(interp_method)(1:mlen)
       else
        write(msg___,667) whoami,methods(interp_method)(1:mlen),whoami,&
     & def
       endif
       do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 666 format(A,' WILL REPARAMETRIZE STRING USING ',A,' INTERPOLATION')
 667 format(A,' WILL REPARAMETRIZE STRING USING ',A,/, &
     &A,' INTERPOLATION TO WITHIN MAX(DS)/MIN(DS) < ',F7.3,' TOLERANCE')
       if (iterations.gt.0) then ; write(msg___,668) whoami, iterations ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 668 format(A,' WITH A MAXIMUM OF ',I5,' ITERATIONS')
       if(interp_method.eq.dst) then ; write(msg___,6680) whoami,dst_cutoff*100.0 ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6680 format(A,' DST INTERPOLATION WILL USE THE LOWER ',F8.4, &
     & '% OF WAVENUMBERS')
      endif
!
! initialize arclength array
      if (.not.allocated(ds)) then
       allocate(ds(nstring-1))
       ds=0.0
      endif
! initialize curvature array
      if (.not.allocated(curv)) then
       allocate(curv(nstring-2))
       curv=0.0
      endif
!
      repa_initialized=.true.
!
      end subroutine smcv_repa_init
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smcv_stat_init(comlyn, comlen)
!
      use sm_var
      use sm_config,only: vtime_offset, rextime_offset
      use cv_common,only:cv_common_neq_work_init, cv_common_rex_read_map
!
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use psf
      use parser
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use mpi
      implicit none
!ccccccccccccccccccccccccccccccccccccccccccccccccc
      CHARACTER(LEN=*) :: COMLYN
      INTEGER :: COMLEN
!
      integer :: klen=0, ierror, i, wtag_len, rex_flen_old
      logical :: found
!
      character(len=80) :: rex_fname_old
!
      character(len=8) :: keyword
      character(len=16) :: whoami
      data whoami/' SMCV_STAT_INIT>'/
!
      logical :: qroot, qprint
 character(len=200) :: msg___(21)=(/'','','','','','','','','','','','','','','','','','','','',''/); integer :: i_
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      qroot=MPI_COMM_STRNG.ne.MPI_COMM_NULL
      qprint=qroot.and.ME_STRNG.eq.0
! begin
! reset iteration counter
! did the user specify it?
      stat_iteration_counter=atoi(get_remove_parameter(comlyn, 'COUN', comlen), -1)
      stat_iteration_counter=max(stat_iteration_counter,0)
      if (stat_iteration_counter.gt.0) then
       if (qprint) then ; write(msg___,639) whoami, stat_iteration_counter ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 639 format(A,' SETTING ITERATION COUNTER TO ',I7)
      endif
!
      cv_fname=''
      output_cv=.false.
      forces_fname=''
      output_forces=.false.
!
      voronoi_fname=''
      output_voronoi_hist=.false.
      output_voronoi_log=.false.
      output_voronoi_map=.false.
!
      rmsd_ave_fname=''
      output_rmsd_ave=.false.
!
      rmsd0_fname=''
      output_rmsd0=.false.
      dsdt_fname=''
      output_dsdt=.false.
!
      c_fname=''
      output_curvature=.false.
!
      s_fname=''
      output_arclength=.false.
!
      fe_fname=''
      output_fe=.false.
!
      work_fname=''
      work_tag=''
      output_work=.false.
!
      rex_fname_old=''
      rex_fname=''
      output_rex_log=.false.
      output_rex_map=.false.
!
      wgt_fname=''
      output_wgt=.false.
!
      M_fname=''
      output_M=.false.
!ccccccccccccccccc first process the RMSD-related commands
!!!!!!!!!!!!!! RMSD from static structure in comp (zts/fts)
      if (remove_tag(comlyn,'RMSD',comlen).gt.0) then ! request for RMSD
       output_rmsd0=.true.
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
       endif
!ccccccccccc print summary
       if (qprint) then
         if (rmsd0_flen.gt.0) then
          write(msg___,660 ) whoami,rmsd0_fname(1:rmsd0_flen)
         else
          write(msg___,661 ) whoami
         endif
         do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
       endif
 660 format(A,' WILL WRITE STRING RMSD TO FILE ',A)
 661 format(A,' WILL WRITE STRING RMSD TO STDOUT.')
!
      endif !! RMSD
!!!!!!!!!!!!!! RMSD from structure at the previous step (zts/fts)
      if (remove_tag(comlyn,'DELS',comlen).gt.0) then
        output_dsdt=.true.
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
        endif
!ccccccccccc print summary
        if (qprint) then
         if (dsdt_flen.gt.0) then
          write(msg___,650 ) whoami,dsdt_fname(1:dsdt_flen)
         else
          write(msg___,651 ) whoami
         endif
         do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
        endif
 650 format(A,' WILL WRITE STRING RMSD(I,I+1) TO FILE ',A)
 651 format(A,' WILL WRITE STRING RMSD(I,I+1) TO STDOUT.')
!
      endif !
!!!!!!!!!!!!!! RMSD from average structure (zts/fts)
      if (remove_tag(comlyn,'RMSA',comlen).gt.0) then
        output_rmsd_ave=.true.
        rmsd_ave_fname=get_remove_parameter(COMLYN,'RANM',COMLEN); rmsd_ave_flen=len_trim(rmsd_ave_fname)
        if (rmsd_ave_flen.eq.0) then
         call warning(whoami, 'NO RMSA FILE NAME SPECIFIED. WILL WRITE TO STDOUT.', 0)
         rmsd_ave_funit=fout
        else
         if (remove_tag(comlyn,'RAAP',comlen).gt.0) then ! APPEND?
           raform='APPEND'
         else
           raform='WRITE'
         endif
        endif
!ccccccccccc print summary
        if (qprint) then
         if (rmsd_ave_flen.gt.0) then
          write(msg___,6500 ) whoami,rmsd_ave_fname(1:rmsd_ave_flen)
         else
          write(msg___,6510 ) whoami
         endif
         do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
        endif
 6500 format(A,' WILL WRITE STRING RMSD FROM AVERAGE STRUC. TO FILE ',A)
 6510 format(A,' WILL WRITE STRING RMSD FROM AVERAGE STRUC. TO STDOUT.')
!
! set number of samples in the average ( to continue a calculation )
        num_average_samples=max(atoi(get_remove_parameter(comlyn, 'NAVE', comlen), 0),0)
! same value for cv
        call cv_set_ave_samples(num_average_samples)
        if (qprint) then ; write(msg___,6511) whoami, num_average_samples ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6511 format(A,' SETTING NUMBER OF SAMPLES IN THE AVERAGE SET TO ',I7)
!
      endif !
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
        endif
!ccccccccccc print summary
        if (qprint) then
         if (s_flen.gt.0) then
          write(msg___,652) whoami,s_fname(1:s_flen)
         else
          write(msg___,653) whoami
         endif
         do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
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
        endif
!ccccccccccc print summary
        if (qprint) then
         if (c_flen.gt.0) then
          write(msg___,6521) whoami,c_fname(1:c_flen)
         else
          write(msg___,6531) whoami
         endif
         do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
        endif
 6521 format(A,' WILL WRITE CURVATURE TO FILE ',A)
 6531 format(A,' WILL WRITE CURVATURE TO STDOUT.')
!
      endif ! CURVATURE
!!!!!!!!!!!!!! FREE ENERGY
      if (remove_tag(comlyn,'FREE',comlen).gt.0) then
        output_fe=.true.
        fe_fname=get_remove_parameter(COMLYN,'FENM',COMLEN); fe_flen=len_trim(fe_fname)
        if (fe_flen.eq.0) then
         call warning(whoami, 'NO F.E. FILE NAME SPECIFIED. WILL WRITE TO STDOUT.', 0)
         fe_funit=fout
        else
         if (remove_tag(comlyn,'FAPP',comlen).gt.0) then ! APPEND?
           feform='APPEND'
         else
           feform='WRITE'
         endif
        endif
!ccccccccccc print summary cccccccccccccccccccccccccccccccccccccc
        if (qprint) then
         if (fe_flen.gt.0) then
          write(msg___,6520) whoami,fe_fname(1:fe_flen)
         else
          write(msg___,6530) whoami
         endif
         do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
        endif
 6520 format(A,' WILL WRITE FREE ENERGY TO FILE ',A)
 6530 format(A,' WILL WRITE FREE ENERGY TO STDOUT.')
!
      endif ! F.E.
!ccccccccccccccccccccccccccc NONEQ. WORK cccccccccccccccccccccccc
      if (remove_tag(comlyn,'WORK',comlen).gt.0) then
        output_work=.true.
        call cv_common_neq_work_init() ! initialize force and position arrays
        work_fname=get_remove_parameter(COMLYN,'WKNM',COMLEN); work_flen=len_trim(work_fname)
        if (work_flen.eq.0) then
         call warning(whoami, 'NO F.E. FILE NAME SPECIFIED. WILL WRITE TO STDOUT.', 0)
         work_funit=fout
        else
         if (remove_tag(comlyn,'WKAP',comlen).gt.0) then ! APPEND?
          wkform='APPEND'
         else
          wkform='WRITE'
         endif
        endif
! specify tag that identifies the work calculated with a particular process
        work_tag=get_remove_parameter(COMLYN,'WTAG',COMLEN); wtag_len=len_trim(work_tag)
        if (wtag_len.eq.0) then
         call warning(whoami, 'WORK TAG NOT SPECIFIED.', 0)
        endif
!ccccccccccc print summary
        if (qprint) then
         if (work_flen.gt.0) then
          write(msg___,6523) whoami,work_fname(1:work_flen)
         else
          write(msg___,6533) whoami
         endif
         do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
        endif
 6523 format(A,' WILL WRITE NON-EQ. WORK TO FILE ',A)
 6533 format(A,' WILL WRITE NON-EQ. WORK TO STDOUT.')
!
      endif ! F.E.
!cccccccccc process CV output options ccccccccccccccccccccccc
      if (remove_tag(comlyn,'COLV',comlen).gt.0) then
! get nergy file name
        cv_fname=get_remove_parameter(COMLYN,'CNAM',COMLEN); cv_flen=len_trim(cv_fname)
!ccccccccccc print summary
        if (cv_flen.gt.0) then
         output_cv=.true.
         if (qprint) then
           write(msg___,6620 ) whoami,cv_fname(1:cv_flen)
           do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
         endif
         if (remove_tag(comlyn,'CVAP',comlen).gt.0) then ! APPEND?
           cvform='APPEND'
         else
           cvform='WRITE'
         endif
        else
          call warning(whoami, 'NO FILE NAME GIVEN. WILL NOT WRITE CV.', 0)
        endif
 6620 format(A,' WILL WRITE CV TIME SERIES TO FILE ',A,'.')
!
      endif ! cv output
!ccccccccccccccccccccccc output weights cccccccccccccccccccccc
      if (remove_tag(comlyn,'WEIG',comlen).gt.0) then
! get nergy file name
        wgt_fname=get_remove_parameter(COMLYN,'WTNM',COMLEN); wgt_flen=len_trim(wgt_fname)
!ccccccccccc print summary
        if (wgt_flen.gt.0) then
         output_wgt=.true.
         if (qprint) then
          write(msg___,6621) whoami,wgt_fname(1:wgt_flen)
          do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
         endif
         if (remove_tag(comlyn,'WTAP',comlen).gt.0) then ! APPEND?
           wgtform='APPEND'
         else
           wgtform='WRITE'
         endif
        else
          call warning(whoami, 'NO FILE NAME GIVEN. WILL NOT WRITE CV WEIGHTS.', 0)
        endif
 6621 format(A,' WILL WRITE CV WEIGHTS TO FILE ',A,'.')
!
      endif ! cv output
!cccccccccccc process Voronoi histogram output options ccccccccccc
      voronoi_flen=0
      if (remove_tag(comlyn,'VORO',comlen).gt.0) then
! get file name
        voronoi_fname=get_remove_parameter(COMLYN,'VNAM',COMLEN); voronoi_flen=len_trim(voronoi_fname)
!ccccccccccc print summary
        if (voronoi_flen.gt.0) then
         output_voronoi_hist=.true.
         if (qprint) then
          write(msg___,6622) whoami,voronoi_fname(1:voronoi_flen)
         endif
        else
         call warning(whoami, 'NO FILE NAME GIVEN. WILL NOT WRITE VORONOI HISTOGRAMS.', 0)
        endif
 6622 format(A,' WILL WRITE VORONOI HISTOGRAMS TO FILE ',A,'.DAT')
!
      endif ! voronoi histograms
!cccccccccccccc voronoi map cccccccccccccccccccccccccccccccccccccccc
      if (remove_tag(comlyn,'VMAP',comlen).gt.0) then
! get file name
        if (voronoi_flen.eq.0) then
         voronoi_fname=get_remove_parameter(COMLYN,'VNAM',COMLEN); voronoi_flen=len_trim(voronoi_fname)
        endif
!
        if (voronoi_flen.gt.0) then
         output_voronoi_map=.true.
         if (qprint) then
          write(msg___,6627) whoami,voronoi_fname(1:voronoi_flen)
          do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
         endif
        else
          call warning(whoami, 'NO FILE NAME GIVEN. WILL NOT WRITE VORONOI MAP.', 0)
        endif
 6627 format(A,' WILL WRITE VORONOI MAP TO FILE ',A,'.MAP')
!
      endif ! voronoi map
!cccccccccccccc voronoi log cccccccccccccccccccccccccccccccccccccccc
      if (remove_tag(comlyn,'VLOG',comlen).gt.0) then
! get file name
        if (voronoi_flen.eq.0) then
         voronoi_fname=get_remove_parameter(COMLYN,'VNAM',COMLEN); voronoi_flen=len_trim(voronoi_fname)
        endif
! check for timestep offset
        vtime_offset=atoi(get_remove_parameter(comlyn, 'VOFF', comlen), 0);
        if (vtime_offset.gt.0) then
         if (qprint) then ; write(msg___,6624) whoami, vtime_offset ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6624 format(A,' WILL OFFSET STEP COUNTER IN VORONOI LOG BY ',I10)
        endif
!
        if (voronoi_flen.gt.0) then
         output_voronoi_log=.true.
         if (qprint) then
          write(msg___,6623) whoami,voronoi_fname(1:voronoi_flen)
          do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
         endif
         if (remove_tag(comlyn,'VLAP',comlen).gt.0) then ! APPEND?
           vlform='APPEND'
         else
           vlform='WRITE'
         endif ! vlap
!
        else
          call warning(whoami, 'NO FILE NAME GIVEN. WILL NOT WRITE VORONOI LOG.', 0)
        endif ! voronoi_flen.gt.0
 6623 format(A,' WILL WRITE VORONOI LOG TO BINARY FILE ',A,'.DAT')
!
      endif ! complete voronoi log
!cccccccccccc replica exchange map cccccccccccccc
      rex_flen=0
      if (remove_tag(comlyn,'REXM',comlen).gt.0) then
! get file name
        rex_fname=get_remove_parameter(COMLYN,'RXNM',COMLEN); rex_flen=len_trim(rex_fname)
! check if user specified an custom map (e.g. from an older run)
        rex_fname_old=get_remove_parameter(COMLYN,'RXOL',COMLEN); rex_flen_old=len_trim(rex_fname_old)
!
        if (rex_flen.gt.0) then
         output_rex_map=.true.
         if (qprint) then
          write(msg___,6721) whoami,rex_fname(1:rex_flen)
          do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
         endif
         if (rex_flen_old.gt.0) then
          if (qprint) then
            write(msg___,6722) whoami,rex_fname_old(1:rex_flen_old) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
            rex_funit=-1
            call files_open(rex_funit, rex_fname_old(1:rex_flen_old), 'FORMATTED', 'READ')
          endif
          call cv_common_rex_read_map(rex_funit)
          if (qprint) then ;call files_close(rex_funit) ; endif
         endif
!
        else
          call warning(whoami, 'NO FILE NAME GIVEN. WILL NOT WRITE REPLICA EXCHANGE MAP.', 0)
        endif
 6721 format(A,' WILL WRITE REPLICA EXCHANGE MAP TO FILE ',A,'.MAP')
 6722 format(A,' WILL RESTART FROM REPLICA EXCHANGE MAP IN FILE ',A)
!
      endif ! replica exchange map
!cccccccccccccc replica exchange log cccccccccccccccccccccccccccccccccccccccc
      if (remove_tag(comlyn,'REXL',comlen).gt.0) then
! get file name
        if (rex_flen.eq.0) then
         rex_fname=get_remove_parameter(COMLYN,'RXNM',COMLEN); rex_flen=len_trim(rex_fname)
        endif
! check for timestep offset
        rextime_offset=atoi(get_remove_parameter(comlyn, 'ROFF', comlen), 0);
        if (rextime_offset.gt.0) then
         if (qprint) then ; write(msg___,6724) whoami, whoami,rextime_offset ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
 6724 format(A,' WILL OFFSET STEP COUNTER IN REPLICA EXCHANGE LOG BY ' &
     & /,A,' ',I10)
        endif
!
        if (rex_flen.gt.0) then
         output_rex_log=.true.
         if (qprint) then
           write(msg___,6723) whoami,whoami,rex_fname(1:rex_flen)
           do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
         endif
         if (remove_tag(comlyn,'RXAP',comlen).gt.0) then ! APPEND?
           rxlform='APPEND'
         else
           rxlform='WRITE'
         endif ! rxap
        else
          call warning(whoami, 'NO FILE NAME GIVEN. WILL NOT WRITE REPLICA EXCHANGE LOG.', 0)
        endif ! rex_flen.gt.0
 6723 format(A,' WILL WRITE REPLICA EXCHANGE LOG TO FILE ',/, &
     & A,' ',A,'.DAT')
!
      endif ! replica exchange log
!cccccccccccccccccc process forces output options cccccccccccccccccc
      if (remove_tag(comlyn,'FORC',comlen).gt.0) then
! get nergy file name
        forces_fname=get_remove_parameter(COMLYN,'FCNM',COMLEN); forces_flen=len_trim(forces_fname)
!ccccccccccc print summary
        if (forces_flen.gt.0) then
         output_forces=.true.
         if (qprint) then
          write(msg___,6625) whoami,forces_fname(1:forces_flen)
         do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
         endif
         if (remove_tag(comlyn,'FCAP',comlen).gt.0) then ! APPEND?
           fform='APPEND'
         else
           fform='WRITE'
         endif
        else
         call warning(whoami, 'NO FILE NAME GIVEN. WILL NOT WRITE AVERAGE FORCE.', 0)
        endif
 6625 format(A,' WILL WRITE AVERAGE FORCE TO FILE ',A,'.')
      endif ! forces
!ccccccccccccccccc process M matrix output options ccccccccccccccccc
      if (remove_tag(comlyn,'MMAT',comlen).gt.0) then
!
        M_fname=get_remove_parameter(COMLYN,'MNAM',COMLEN); M_flen=len_trim(M_fname)
        if (M_flen.gt.0) then
         output_M=.true.
         if (qprint) then ; write(msg___,6626) whoami,M_fname(1:M_flen) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ; endif
        else
         call warning(whoami, 'NO FILE NAME GIVEN. WILL NOT WRITE M TENSOR.', 0)
        endif
 6626 format(A,' WILL WRITE TENSOR M TO FILE ',A,'.')
      endif
!
! if we got this far, we are probably OK
      stat_initialized=.true.
!
      end subroutine smcv_stat_init
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smcv_stat()
! output statistics for SMCV
      use sm_var
      use sm_config, only: calc_cv_para, calc_Mtensor_para
      use cv_common
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use parser
      use mpi ! deal with other platforms later (never)
!
      implicit none
!
!
      integer :: me, ierror, i, ifile, k, fmt_r_len
      real*8 :: rmsd0, rmsd0_all(nstring), dsdt, dsdt_all(nstring)
      real*8 :: mywork, allwork(nstring)
      character(len=8) :: dummy, work_tags(nstring)
      character(len=80) :: fmt, fmt_real, fmt_int ! format strings for output
      integer :: oldiol
      character(len=11) :: whoami
      data whoami/' SMCV_STAT>'/
!
      logical :: qroot, qprint, qgrp
!
      qroot=MPI_COMM_STRNG.ne.MPI_COMM_NULL
      qprint=qroot.and.ME_STRNG.eq.0
      qgrp=MPI_COMM_LOCAL.ne.MPI_COMM_NULL.and.SIZE_LOCAL.gt.1
!
! ad hoc fix for REX
!cccccccccccccccccccccc begin ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! check if the user has made an initialization call
      if (.not.smcv_initialized) call smcv_init()
      if (.not.stat_initialized) then
       call warning(whoami, 'NO OUTPUT OPTIONS SELECTED. NOTHING DONE', 0)
       return
      endif
!
      stat_iteration_counter=stat_iteration_counter+1
! define number format string for output
!
      write(fmt_real,*) nstring
      fmt_r_len=len(fmt_real)
      fmt_r_len=min(max(0,fmt_r_len),len(fmt_real));fmt_real(fmt_r_len+1:)='';call adjustleft(fmt_real,(/' ',tab/));fmt_r_len=len_trim(fmt_real)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_rmsd0) then
! weighting is taken care of in cv_common class
          rmsd0=cv_common_rmsd(1,3) ! we will keep the initial reference coords in col. 3
          if (qroot) call mpi_gather(rmsd0,1,MPI_REAL & ! heads communicate
     & ,rmsd0_all,1,MPI_REAL,0, &
     & MPI_COMM_STRNG, ierror)
          if (qprint) then ! root writes
           if (rmsd0_funit.eq.fout) then
            fmt='("RMSD0> ",I5," ",'//fmt_real(1:fmt_r_len)//'F15.5)'
           else
            rmsd0_funit=-1
            call files_open(rmsd0_funit, rmsd0_fname, 'FORMATTED', rform)
            fmt='(I5," ",'//fmt_real(1:fmt_r_len)//'F15.5)'
           endif
           write(rmsd0_funit,fmt) stat_iteration_counter, &
     & (rmsd0_all(i),i=1,nstring)
!
           if (rmsd0_funit.ne.fout) then
            call files_close(rmsd0_funit)
           endif
          endif ! qprint
          rform='APPEND'
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_dsdt) then
! weighting is taken care of in cv_common class
          dsdt=cv_common_rmsd(1,2) ! we will keep the previous coords in col. 2
          if (qroot) call mpi_gather(dsdt,1,MPI_REAL &
     & ,dsdt_all,1,MPI_REAL,0, &
     & MPI_COMM_STRNG, ierror)
          if (qprint) then ! root writes
           if (dsdt_funit.eq.fout) then
            fmt='("DLEN> ",I5," ",'//fmt_real(1:fmt_r_len)//'F15.5)'
           else
            dsdt_funit=-1
            call files_open(dsdt_funit, dsdt_fname, 'FORMATTED', dform)
            fmt='(I5," ",'//fmt_real(1:fmt_r_len)//'F15.5)'
           endif
           write(dsdt_funit,fmt) stat_iteration_counter, &
     & (dsdt_all(i),i=1,nstring)
! flush unit: close and reopen
           if (dsdt_funit.ne.fout) then
            call files_close(dsdt_funit)
           endif
          endif ! qprint
          dform='APPEND'
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_rmsd_ave) then ! rmsd with respect to the running average structure
! weighting is taken care of in cv_common class
! call to update running average
          call cv_common_update_ave()
          rmsd0=cv_common_rmsd(1,4) ! we will keep the average coords in col. 4
          if (qroot) call mpi_gather(rmsd0,1,MPI_REAL &
     & ,rmsd0_all,1,MPI_REAL,0, &
     & MPI_COMM_STRNG, ierror)
          if (qprint) then ! root writes
           if (rmsd_ave_funit.eq.fout) then
            fmt='("RMSD_AVE> ",I5," ",'//fmt_real(1:fmt_r_len)//'F15.5)'
           else
            rmsd_ave_funit=-1
            call files_open(rmsd_ave_funit, rmsd_ave_fname, 'FORMATTED', raform)
            fmt='(I5," ",'//fmt_real(1:fmt_r_len)//'F15.5)'
           endif
           write(rmsd_ave_funit,fmt) stat_iteration_counter, &
     & (rmsd0_all(i),i=1,nstring)
! flush unit: close and reopen
           if (rmsd_ave_funit.ne.fout) then
            call files_close(rmsd_ave_funit)
           endif
          endif ! qprint
          raform='APPEND'
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_arclength) then
       if (qprint) then
         write(fmt_int,'(I5)') stat_iteration_counter
         if (s_funit.eq.fout) then
          fmt='("ARCL> '//fmt_int(1:5)//' ",'                           &
     & //fmt_real(1:fmt_r_len)//'F15.5)'
         else
          s_funit=-1
          call files_open(s_funit, s_fname, 'FORMATTED', sform)
          fmt='("'//fmt_int(1:5)//' ",'//fmt_real(1:fmt_r_len)//'F15.5)'
         endif
         call cv_common_print_ds(s_funit, fmt)
! flush unit: close and reopen
         if (s_funit.ne.fout) then
          call files_close(s_funit)
         endif
! done
       endif ! qprint
       sform='APPEND'
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_curvature) then
       if (qprint) then
         write(fmt_int,'(I5)') stat_iteration_counter
         if (c_funit.eq.fout) then
          fmt='("CURV> '//fmt_int(1:5)//' ",'                           &
     & //fmt_real(1:fmt_r_len)//'F11.5)'
         else
          c_funit=-1
          call files_open(c_funit, c_fname, 'FORMATTED', cform)
          fmt='("'//fmt_int(1:5)//' ",'//fmt_real(1:fmt_r_len)//'F11.5)'
         endif
         call cv_common_print_curvature(c_funit, fmt)
! flush unit: close and reopen
         if (c_funit.ne.fout) then
          call files_close(c_funit)
! done
         endif
       endif ! qprint
       cform='APPEND'
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_fe) then
       if (qprint) then
         write(fmt_int,'(I5)') stat_iteration_counter
         if (fe_funit.eq.fout) then
          fmt='("FE> '//fmt_int(1:5)//' ",'                             &
     & //fmt_real(1:fmt_r_len)//'F15.5)'
         else
          fe_funit=-1
          call files_open(fe_funit, fe_fname, 'FORMATTED', feform)
          fmt='("'//fmt_int(1:5)//' ",'//fmt_real(1:fmt_r_len)//'F15.5)'
         endif
         call cv_common_print_feav(fe_funit, fmt)
! flush unit: close and reopen
         if (fe_funit.ne.fout) then
          call files_close(fe_funit)
! done
         endif
       endif ! qprint
       feform='APPEND'
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_work) then
! gather work
       mywork=cv_common_neq_get_work()
! if running in parallel, need to reduce work from slave nodes
       if (MPI_COMM_LOCAL.ne.MPI_COMM_NULL.and.SIZE_LOCAL.gt.1 &
     & .and.calc_cv_para) then
        call MPI_REDUCE(mywork, allwork(1), 1, MPI_REAL, &
     & MPI_SUM, 0, MPI_COMM_LOCAL, ierror) ! reduce on all group roots
        mywork=allwork(1)
       endif
! gather work from all nodes into one output buffer
       if (qroot) then
        call mpi_gather(mywork, 1, MPI_REAL, &
     & allwork, 1, MPI_REAL, &
     & 0, MPI_COMM_STRNG, ierror)
        call mpi_gather(work_tag, 8, MPI_BYTE, &
     & work_tags, 8, MPI_BYTE, 0, MPI_COMM_STRNG, &
     & ierror)
       endif
! write
       if (qprint) then
         write(fmt_int,'(I5)') stat_iteration_counter
         if (work_funit.eq.fout) then
          fmt='("WORK> '//fmt_int(1:5)//' ",A8,F15.5)'
         else
          work_funit=-1
          call files_open(work_funit, work_fname, 'FORMATTED', wkform)
          fmt='("'//fmt_int(1:5)//' ",A8,F15.5)'
         endif
         do i=1,nstring
          write(work_funit,fmt) work_tags(i), allwork(i)
         enddo
! flush unit: close and reopen
         if (work_funit.ne.fout) then
          call files_close(work_funit)
         endif
! done
       endif ! qprint
       wkform='APPEND'
      endif ! output_work
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_voronoi_hist) then ! output voronoi data
        if (voronoi_flen.eq.0) then
         call warning(whoami, 'NO FILE NAME SPECIFIED. WILL NOT WRITE VORONOI DATA.', 0)
        else
         if (qprint) then
          ifile=-1
          voronoi_fname(voronoi_flen+1:voronoi_flen+4)='.dat'
          call files_open(ifile, voronoi_fname(1:voronoi_flen+4), 'FORMATTED', 'WRITE')
          voronoi_fname(voronoi_flen+1:)=''
         endif
         call cv_common_print_voro_data(ifile) ! all root processes enter
         if (qprint) then ; call files_close(ifile) ; endif
        endif
      endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_voronoi_log) then
       if (voronoi_flen.eq.0) then
         call warning(whoami, 'NO FILE NAME SPECIFIED. WILL NOT WRITE VORONOI LOG.', 0)
       else
         if (qprint) then
           vlog_funit=-1
           voronoi_fname(voronoi_flen+1:voronoi_flen+4)='.log'
           call files_open(vlog_funit, voronoi_fname(1:voronoi_flen+4), 'UNFORMATTED', vlform)
           voronoi_fname(voronoi_flen+1:)=''
         endif
         vlform='APPEND'
         if (qroot) call cv_common_voronoi_print_log(vlog_funit)
! flush unit:
         if (qprint) then ; call files_close(vlog_funit) ; endif
       endif
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_voronoi_map) then ! output voronoi map
        if (voronoi_flen.eq.0) then
         call warning(whoami, 'NO FILE NAME SPECIFIED. WILL NOT WRITE VORONOI MAP.', 0)
        else
! put 'whereami' into the map
          if (qroot.and.SIZE_STRNG.gt.1) then
           call MPI_ALLGATHER(cv%voronoi_whereami, 1, MPI_INTEGER, &
     & cv%voronoi_map, 1, MPI_INTEGER, MPI_COMM_STRNG, ierror)
          else
           cv%voronoi_map(mestring+1)=cv%voronoi_whereami
          endif
          if (qgrp) then
           call mpi_bcast(cv%voronoi_map,nstring,MPI_INTEGER,0,MPI_COMM_LOCAL,ierror)
          endif ! qgrp
!
         if (qroot) then
          ifile=-1
          voronoi_fname(voronoi_flen+1:voronoi_flen+4)='.map'
          if (qprint) then
           call files_open(ifile, voronoi_fname(1:voronoi_flen+4), 'FORMATTED', 'WRITE')
           voronoi_fname(voronoi_flen+1:)=''
          endif
          call cv_common_print_voro_map(ifile)
          if (qprint) then ; call files_close(ifile) ; endif
         endif ! qroot
        endif
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_rex_map) then ! output replica exchange map
        if (rex_flen.eq.0) then
         call warning(whoami, 'NO FILE NAME SPECIFIED. WILL NOT WRITE REPLICA EXCHANGE MAP.', 0)
        else
         if (qprint) then
          rex_funit=-1
          rex_fname(rex_flen+1:rex_flen+4)='.map'
          call files_open(rex_funit, rex_fname(1:rex_flen+4), 'FORMATTED', 'WRITE')
          rex_fname(rex_flen+1:)=''
!
          call cv_common_rex_print_map(rex_funit)
!
          call files_close(rex_funit)
         endif ! qprint
        endif
      endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_rex_log) then
       if (rex_flen.eq.0) then
         call warning(whoami, 'NO FILE NAME SPECIFIED. WILL NOT WRITE REPLICA EXCHANGE LOG.', 0)
       else
        if (qprint) then
         rex_funit=-1
         rex_fname(rex_flen+1:rex_flen+4)='.dat' ! append to name
         call files_open(rex_funit, rex_fname(1:rex_flen+4), 'FORMATTED', rxlform)
         rex_fname(rex_flen+1:)='' ! erase extension
        endif
        rxlform='APPEND'
!
        call cv_common_rex_print_log(rex_funit)
! flush unit:
        if (qprint) then ; call files_close(rex_funit) ; endif
       endif
      endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_cv) then
        if (qprint) then
         cv_funit=-1
         call files_open(cv_funit, cv_fname, 'FORMATTED', cvform)
         write(cv_funit,'("% ",I8)') stat_iteration_counter ! % is a MATLAB comment
        endif
        cvform='APPEND'
        call cv_common_print_global(cv_funit)
! flush unit:
        if (qprint) then ; call files_close(cv_funit) ; endif
      endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_forces) then
       if (qprint) then
        forces_funit=-1
        call files_open(forces_funit, forces_fname, 'FORMATTED', fform);
        write(forces_funit,'("% ",I8)') stat_iteration_counter ! % is a MATLAB comment
       endif
       fform='APPEND'
       call cv_common_print_forces(forces_funit)
! flush unit:
       if (qprint) then ; call files_close(forces_funit) ; endif
      endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_wgt) then
       if (qprint) then
         wgt_funit=-1
         call files_open(wgt_funit, wgt_fname, 'FORMATTED', wgtform);
         write(wgt_funit,'("% ",I8)') stat_iteration_counter ! % is a MATLAB comment
         call cv_common_print_wgt(wgt_funit)
! flush unit:
         call files_close(wgt_funit)
       endif
       wgtform='APPEND'
      endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (output_M) then
! if running in parallel, combine partial M entries
       if (qgrp.and.calc_Mtensor_para) then
! call MPI_ALLREDUCE(cv%M(1,1,1),cv%M(1,1,2),max_cv_common*cv%num_cv, ! will broadcast all rows, but only num_cv columns
! & MPI_REAL, MPI_SUM, MPI_COMM_LOCAL, ierror)
       call MPI_REDUCE(cv%M(1,1,1),cv%M(1,1,2),max_cv_common*cv%num_cv, & ! will broadcast all rows, but only num_cv columns
     & MPI_REAL, MPI_SUM, 0, MPI_COMM_LOCAL, ierror)
       else ! qgrp
         cv%M(1:cv%num_cv,1:cv%num_cv,2)=cv%M(1:cv%num_cv,1:cv%num_cv,1)
       endif ! qgrp
!
       if (qroot) then
        if (M_flen.eq.0) then
         call warning(whoami, 'NO FILE NAME SPECIFIED. WILL NOT WRITE M TENSOR.', 0)
        else
         ifile=-1
         call files_open(ifile, M_fname(1:M_flen), 'FORMATTED', 'WRITE')
         call cv_common_print_M_global(ifile)
         call files_close(ifile)
        endif ! M_flen
       endif ! qroot
      endif ! output_M
! ad hoc fix for REX
!
      end subroutine smcv_stat
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smcv_repa()
! this is routine is just a wrapper, so that any subroutine can call a "default" repa.
      use sm_var,only: interp_method, iterations, def, dst_cutoff, &
     & repa_initialized
      use cv_common,only: cv_common_repa
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use parser
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use mpi
!
      implicit none
 character(len=200) :: msg___(21)=(/'','','','','','','','','','','','','','','','','','','','',''/); integer :: i_
! local variables
      character(len=11) :: whoami
      logical :: qprint
      data whoami/' SMCV_REPA>'/
! check if the user has made an initialization call
!
      qprint=(MPI_COMM_STRNG.ne.MPI_COMM_NULL.and.ME_STRNG.eq.0)
!
      if (.not.repa_initialized) then
       call warning(whoami, 'NO REPARAMETRIZATION OPTIONS SELECTED. NOTHING DONE.', 0)
       return
      endif
      if (qprint) then ; write(msg___,690) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_),5);enddo;msg___='' ; endif
 690 format(/A,' CALLING STRING REPARAMETRIZATION.')
      call cv_common_repa(interp_method,def,iterations,dst_cutoff)
!
      end subroutine smcv_repa
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cv_set_ave_samples(n)
      use cv_common,only: cv_common_set_ave_samples
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      implicit none
      integer :: n
      call cv_common_set_ave_samples(n)
      end subroutine cv_set_ave_samples
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!**CHARMM_ONLY**!##ENDIF
