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
module output
 implicit none
 integer, save :: fout=-1
 logical, save :: output_initialized=.false.
 logical, save :: qprint = .true.
 logical, private, save :: loud = .true.
 character(LEN=80) :: msg__ ! buffer for output (see message.src) macro
 character(LEN=10), public, parameter :: realfmt = 'ES20.7E3'
!
 integer, private, parameter :: minerrorlev=0, minwarnlev=0
 integer, private :: l
 integer, private, save :: warnlev, errorlev ! last warning/error levels
 character*9, private :: stat
 contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine error(whoami,msg,level)
  integer :: level
  character*(*) :: whoami, msg
!
  if (.not.output_initialized) call output_init()
  if (level.lt.minerrorlev) then
   stat=' FATAL'; l=6
  else
   stat=' NONFATAL'; l=9
  endif
!
  if (qprint.and.loud) then
   write(fout,'(5A)') stat(1:l),' ERROR (',whoami,'): ', msg
  endif
!
  if (level.lt.minerrorlev) stop
!
 end subroutine error
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine terminate(whoami)
  character*(*) :: whoami
!
  if (.not.output_initialized) call output_init()
!
  if (qprint) then
   write(fout,'(5A)') ' TERMINATION INVOKED FROM (',whoami,').'
  endif
!
  stop
!
 end subroutine terminate
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine warning(whoami,msg,level)
  integer :: level
  character*(*) :: whoami, msg
!
  if (.not.output_initialized) call output_init()
  if (level.lt.minwarnlev) then
   stat=' FATAL'; l=6
  else
   stat=' NONFATAL'; l=9
  endif
!
  if (qprint) then
   write(fout,'(5A)') stat(1:l),' WARNING (',whoami,'): ', msg
  endif
!
  warnlev=min(warnlev,level) ! keep the most severe level
!
! termination due to a fatal warning is triggered by user (usually via routine below)
!
 end subroutine warning
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function fatal_warning(comm)
 integer*4 :: comm ! communicator
 integer :: level, bug
 logical :: fatal_warning
!
!
! compute the maximum value of errorcodes across all nodes; then decice whether to terminate
!
 level=warnlev ! current value of the error level
!
!
 fatal_warning=(level.lt.minerrorlev)
!
 end function fatal_warning
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine message(whoami,msg)
  character*(*) :: whoami, msg
!
  if (.not.output_initialized) call output_init()
!
  if (qprint) then
   write(fout,'(4A)') ' MESSAGE (',whoami,'): ', msg
  endif
!
 end subroutine message
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine plainmessage(msg)
! write message without ay other info.
  character*(*) :: msg
!
  if (.not.output_initialized) call output_init()
!
  if (qprint) then
   write(fout,'(4A)') msg
  endif
!
 end subroutine plainmessage
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
 subroutine output_init(filename)
 character*(*), optional :: filename
 character*400 :: fname
 integer :: flen
 call output_done()
!
 if (present(filename)) then
  fname=adjustl(filename)
  flen=len_trim(fname)
  fout=987
  if (qprint) &
& open(unit=fout,file=fname(1:flen),form='FORMATTED',status='UNKNOWN')
 else
  fout=6
 endif
!
 output_initialized=.true.
 end subroutine output_init
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
 subroutine output_done()
 if (output_initialized) then
  if (qprint.and.fout.ne.0.and.fout.ne.6.and.fout.ne.5) close(fout)
  output_initialized=.false.
  fout=-1
 endif
 end subroutine output_done
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module output
