/*COORDINATES AND MASSES:*/
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
module files
 use output, only: message, warning
 use ivector
 use parser, only: toupper, itoa
 implicit none
 public files_open
 public files_close
 private files_initialize
 public files_done
 private files_assign_handle
 type (int_vector), private, save :: handles
 logical, public, save :: files_initialized=.false.
 integer, private :: ihandle ! initial file handle for automatic assignment
 integer, parameter, public :: reserved_handles(3)=(/0, 5, 6/) ! disallowed handles
 integer, parameter, public :: ihandle_max=999 ! maximum unit number
 integer, parameter, public :: ihandle_min=101 ! minimum unit number
 contains
!******************************************************************
 subroutine files_open(handle_, name_, form_, action_)
 integer :: handle_
 character(len=*), intent(in) :: name_, form_, action_
 character(len=len(form_)) :: form__
 character(len=len(action_)) :: action__
 character(len=132) :: name__
 character(len=7) :: msg__, stat_
 logical :: exist_, opened_, qread, qwrite, qappend
 character*10, parameter :: whoami='FILES_OPEN'
 integer :: ahandle, ioerr, i
!
 form__=form_
 action__=action_
! access checks
 call toupper(form__)
 call toupper(action__)
 qread=.false.; qwrite=.false.; qappend=.false.
!
 if (action__.eq.'READ') then ; qread=.true. ; msg__='reading' ; stat_='OLD'
 elseif (action__.eq.'WRITE') then ; qwrite=.true.; msg__='writing'; stat_='UNKNOWN'
 elseif (action__.eq.'APPEND') then ; qappend=.true.; msg__='appending'; stat_='OLD'
 endif
!
! write(0,*) handle_
! write(0,*) name_
! write(0,*) form_
!
 if (.not.files_initialized) call files_initialize()
! check that file exists
 inquire(file=name_,exist=exist_, opened=opened_, number=ahandle)
 if (.not.exist_.and.(qread.or.qappend)) then
  call warning(whoami,'Cannot find file '//name_(1:len_trim(name_)),0)
  handle_=-1
 else ! file exists, or we are writing
  if (opened_) then
   call warning(whoami,'File '//name_(1:len_trim(name_))//' is already connected to unit '//itoa(ahandle)//'. Disconnecting.',0)
   call files_close(ahandle)
  endif
!
  if (any(reserved_handles.eq.handle_)) then ; call warning(whoami, 'Invalid file handle requested',0); handle_=-1
  elseif(.not. (form__.eq.'FORMATTED'.or.form__.eq.'UNFORMATTED')) then
   call warning(whoami,'File format '//form__(1:len_trim(form__))//' is not allowed.',0); handle_=-1
  elseif(.not.(qread.or.qwrite.or.qappend)) then
   call warning(whoami,'File access '//action__(1:len_trim(action__))//' is not allowed.',0); handle_=-1
  elseif (handle_.lt.0) then ; handle_=files_assign_handle()
  else ! default: check
   inquire(handle_, name=name__, opened=opened_)
   if (opened_) then
    call warning(whoami,'Unit '//itoa(handle_)//' is already connected to file '//name__(1:len_trim(name__))//'. Disconnecting.',0)
    call files_close(handle_)
   endif
  endif
 endif ! exist_
!
 if (handle_.gt.-1) then
! open file
   call message(whoami, 'Opening '//form__(1:len_trim(form__))//' file "'//name_(1:len_trim(name_))//'" for '//msg__//'.')
!
   if (qappend) then
    open(unit=handle_, file=name_, form=form__, status=stat_, iostat=ioerr, action='WRITE', position=action__)
   else
    open(unit=handle_, file=name_, form=form__, status=stat_, iostat=ioerr, action=action__)
   endif
!
   if (ioerr.ne.0) then
    call warning(whoami, 'File open failed with error code '//itoa(ioerr),0)
    handle_=-1
   else ! all OK
    i=int_vector_uadd(handles, handle_)
   endif
!
 endif ! handle
!
 end subroutine files_open
!******************************************************************
 subroutine files_initialize()
 character(len=16), parameter :: whoami='FILES_INITIALIZE'
 if (files_initialized) call files_done()
 call int_vector_init(handles)
 if (.not.handles%initialized) then
  call warning(whoami,'Cannot initialize vector of file handles',-1)
  return
 else
  files_initialized=.true.
  ihandle=ihandle_min
 endif
 end subroutine files_initialize
!******************************************************************
 subroutine files_close(handle_)
 integer, intent(inout) :: handle_
 integer :: i
 logical :: opened_, ok
 character(len=11), parameter :: whoami='FILES_CLOSE'
!
 if (.not.files_initialized) call files_initialize()
!
 i=int_vector_getind(handles, handle_)
 if (i.gt.0.and..not.any(reserved_handles.eq.i)) then ! valid entry; check using inquire next
  inquire(handle_, opened=opened_)
  if (.not.opened_) then
   call warning(whoami,'Unit number is valid, but is not open (according to INQUIRE)',0)
  else
   call message(whoami,'Closing unit '//itoa(handle_))
   close(handle_)
   handle_=-1
  endif
  ok=int_vector_delete(handles, i)
  if (.not.ok) call warning(whoami,'Cannot close unit '//itoa(handle_)//' (internal error).',-1)
 endif
!
 end subroutine files_close
!******************************************************************
 subroutine files_done()
 integer :: i
 if (handles%initialized) then
  do i=1, handles%last
   call files_close(handles%i(i))
  enddo
  call int_vector_done(handles)
 endif
 files_initialized=.false.
!
 end subroutine files_done
!******************************************************************
 function files_assign_handle()
 integer :: files_assign_handle
 logical :: opened_
 character(len=19), parameter :: whoami='FILES_ASSIGN_HANDLE'
!
 if (.not.files_initialized) call files_initialize()
 do
  if (ihandle.ge.ihandle_max) then
  call warning(whoami,'Maximum number of allowed opened units ('//itoa(ihandle_max-ihandle_min+1)//') exceeded',-1);ihandle=-1;exit
  else
   inquire(ihandle, opened=opened_)
   if (.not.opened_) exit
  endif
  ihandle=ihandle+1
 enddo
 files_assign_handle=ihandle
!
 end function files_assign_handle
!******************************************************************
end module files
