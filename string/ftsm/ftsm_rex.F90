! **********************************************************************!
! This source file was was generated automatically from a master source !
! code tree, which may or may not be distributed with this code, !
! because it is up to the distributor, and not up to me. !
! If you edit this file (rather than the master source file) !
! your changes will be lost if another pull from the master tree occurs.!
! In case you are wondering why, this approach makes it possible for !
! me to have the same master source code interfaced with different !
! applications (some of which are written in a way that is quite far !
! from being object-oriented) at the source level. !
! **********************************************************************!
!
! FTSM_REX.MOD
!
! REPLICA EXCHANGE MODULE FOR THE FINITE TEMPERATURE STRING METHOD
!**CHARMM_ONLY**!##IF STRINGM
!
      module ftsm_rex
      use ivector
      use ftsm_var, only: nstring, ftsm_initialized
!
     
!
      private
!
      real*8, save, public :: rex_beta
      integer, save, allocatable, public :: rex_map(:)
      type (int_vector), save, public :: rex_log
!
      logical, save, public :: rex_initialized=.false.
!
      public ftsm_rex_init
      public ftsm_rex_done
      public ftsm_rex_set_temp
      public ftsm_rex_print_map
      public ftsm_rex_read_map
      public ftsm_rex_print_log
!
      contains
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine ftsm_rex_init(temp)
      use constants
!
       real*8, optional :: temp
       real*8 :: t
       integer :: i
!
       if (.not.rex_initialized) then
        if (present(temp)) then ; t=temp; else ; t=300d0; endif
        if (t.gt.0) rex_beta=1d0/(t*kboltzmann)
        if (.not.ftsm_initialized) return
!
        allocate(rex_map(nstring))
        rex_map=(/ (i, i=0,nstring-1) /)
        call int_vector_init(rex_log)
        rex_initialized=.true.
       endif
!
       end subroutine ftsm_rex_init
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine ftsm_rex_done()
       if (rex_initialized) then
        deallocate(rex_map)
        call int_vector_done(rex_log)
        rex_initialized=.false.
       endif
       end subroutine ftsm_rex_done
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine ftsm_rex_set_temp(temp)
      use constants
       real*8, optional :: temp
       real*8 :: t
!
       if (.not.rex_initialized) call ftsm_rex_init()
       if (present(temp)) then ; t=temp; else ; t=300d0; endif
       if (t.gt.0) rex_beta=1d0/(t*kboltzmann)
       end subroutine ftsm_rex_set_temp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine ftsm_rex_print_map(iunit,fmt)
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning
! only root process should call
       integer :: iunit
       character(len=*), optional :: fmt
! local
       integer :: i
       character(80) :: frm
       character(len=20) :: whoami
       data whoami /' FTSM_PRINT_REX_MAP>'/
! begin
       if (.not.rex_initialized) then
        write(0,*) 'WARNING FROM: ',whoami,': ','REX NOT INITIALIZED. NOTHING DONE.'
        return
       endif
!
       if (.not.present(fmt)) then
        write(frm,'("(",I5,"I5)")') nstring
       else
        frm=fmt
       endif
       write(iunit,frm) (/ (i, i=0,nstring-1) /)
       write(iunit,frm) rex_map(1:nstring)
       end subroutine ftsm_rex_print_map
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine ftsm_rex_read_map(iunit)
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning
      use multicom_aux
      use mpi
       integer :: iunit, ierror
       character(len=19) :: whoami
       data whoami /' FTSM_READ_REX_MAP>'/
! begin
       if (.not.rex_initialized) then
        call ftsm_rex_init()
       endif
!
       if (MPI_COMM_STRNG.ne.MPI_COMM_NULL) then
        if (ME_STRNG.eq.0) then
         read(iunit,*) rex_map(1:nstring) ! first row contains indices 0 -- nstring-1
         read(iunit,*) rex_map(1:nstring) ! second row is what we want
         if (any(rex_map.lt.0)) write(0,*) 'WARNING FROM: ',whoami,': ','READ NEGATIVE RANK.'
        endif ! ME_
        if (SIZE_STRNG.gt.1) call mpi_bcast(rex_map,nstring,MPI_INTEGER,0,MPI_COMM_STRNG,ierror)
       endif ! MPI_COMM
! broadcast to slave nodes
       if (ME_LOCAL.ne.MPI_UNDEFINED.and.SIZE_LOCAL.gt.1) &
! & call MPI_BCAST(cv%rex_map, nstring, MPI_INTEGER,
! & 0,MPI_COMM_LOCAL,ierr)
     & call mpi_bcast(rex_map,nstring,MPI_INTEGER,0,MPI_COMM_LOCAL,ierror)
!
       end subroutine ftsm_rex_read_map
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine ftsm_rex_print_log(iunit, fmt)
! assume that unit is prepared
! NOTE that this is a global print!
      use multicom_aux
      use mpi
!
       integer :: iunit
       character(len=*), optional :: fmt
! local
       character(80) :: frm
       integer :: i
       integer*4 :: rex_log_size4(nstring)
       integer :: rex_log_size8(nstring)
       integer*4 :: rex_log_disp4(nstring)
       integer :: total_size
       integer, allocatable, dimension(:) :: rex_log_all
       integer :: ierror, type
       logical :: qroot
!
! do work
! gather all data on root
       qroot=MPI_COMM_STRNG.ne.MPI_COMM_NULL
!
       if (.not.rex_initialized) call ftsm_rex_init()
!
       if (qroot.and.SIZE_STRNG.gt.1) then
! calculate size of logs
        rex_log_size8=0
! call MPI_ALLGATHER(rex_log%last,1,type,
! & rex_log_size8,1,type,
! & MPI_COMM_STRNG,error)
        call MPI_GATHER(rex_log%last,1,MPI_INTEGER, &
     & rex_log_size8,1,MPI_INTEGER, &
     & 0,MPI_COMM_STRNG,ierror)
        call mpi_bcast(rex_log_size8,nstring,MPI_INTEGER,0,MPI_COMM_STRNG,ierror)
!
        total_size=sum(rex_log_size8)
        rex_log_size4=rex_log_size8 ! type cast to 4 byte integer
! allocate space to hold entire log
        allocate(rex_log_all(total_size))
! calculate send displacements
        rex_log_disp4(1)=0;
        do i=1,SIZE_STRNG-1
         rex_log_disp4(i+1)=rex_log_disp4(i)+rex_log_size4(i)
        enddo
! now gather the logs
! call MPI_ALLGATHERV(rex_log%i,rex_log%last,type,
! & rex_log_all,rex_log_size4,rex_log_disp4,type,
! & MPI_COMM_STRNG,ierror)
        call MPI_GATHERV(rex_log%i,rex_log%last,type, &
     & rex_log_all,rex_log_size4,rex_log_disp4,type, &
     & 0,MPI_COMM_STRNG,ierror)
!
        if (.not.present(fmt)) then
         frm='(2I5,I8)'
        else
         frm=fmt
        endif
!
      if(ME_STRNG.eq.0.and.total_size.gt.0) write(iunit,frm) rex_log_all
!
        call int_vector_init(rex_log) ! erase log
        deallocate(rex_log_all)
       endif ! STRNG
       end subroutine ftsm_rex_print_log
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end module ftsm_rex
!
!**CHARMM_ONLY**!##ENDIF
