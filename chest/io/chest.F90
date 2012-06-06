/*#define __WRN(__WHO,__MSG) write(0,*) 'WARNING FROM: ',__WHO,': ',__MSG*/
/*#define __PRINT(__MSG) write(0,'(A)') __MSG*/
/*COORDINATES AND MASSES:*/
/*#define __INDX(__STR, __STRLEN, __TEST, __TESTLEN)  index(__STR(1:min(__STRLEN,len(__STR))),__TEST(1:min(__TESTLEN,len(__TEST))))*/
/*#define __WRN(__WHO,__MSG) write(0,*) 'WARNING FROM: ',__WHO,': ',__MSG*/
/*#define __PRINT(__MSG) write(0,'(A)') __MSG*/
/*COORDINATES AND MASSES:*/
/*#define __INDX(__STR, __STRLEN, __TEST, __TESTLEN)  index(__STR(1:min(__STRLEN,len(__STR))),__TEST(1:min(__TESTLEN,len(__TEST))))*/
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine chest_read_scalar(filename,q,nx,ny,nz,qbin) ! scalar function
  use output
  use files
  use gridsize, only : me, communicator
  use parser, only: adjustleft
  implicit none
!
  real*8 :: q(nx,ny,nz)
  integer :: nx, ny, nz
  integer :: mx, my, mz
  integer :: nblock
  logical :: qbin
  character(len=11) :: fmt
  character(len=17) , parameter :: whoami='CHEST_READ_SCALAR'
!
  logical :: bin
  logical :: q3d
!
  integer :: i, j, k, ioerr
!
  character(len=*) :: filename
  character(len=80) :: fname
  integer :: flen, fid=-1
!
  bin=qbin
  if (bin) then; fmt='UNFORMATTED'; else ; fmt='FORMATTED'; endif
!
  fname=filename
  call adjustleft(fname)
  flen=len_trim(fname)
  if (flen.gt.0) then
   if (me.eq.0) then
    call files_open(fid, name_=fname(1:flen), form_=fmt, action_='READ')
    if (fid.lt.0) then
     call warning(whoami, 'Cannot open input file. Abort.',-1)
    endif
   endif ! me
   if (fatal_warning()) return
  else
   call warning(whoami, 'Input file name not specified. Abort.',0)
   return
  endif
!
  q3d=(nz.gt.3)
  mx=nx; my=ny; mz=nz ! initialize
!
  if (me.eq.0) then ! may change parallel behavior in the future
!************** header *******
!
  if (bin) then
   if (q3d) then
    read(fid, IOSTAT=ioerr) mx, my, mz
   else
    read(fid, IOSTAT=ioerr) mx, my
   endif ! q3d
  else ! qbin
   if (q3d) then
    read(fid,*, IOSTAT=ioerr) mx, my, mz
   else
    read(fid,*, IOSTAT=ioerr) mx, my
   endif ! q3d
  endif ! qbin
!
  if (ioerr.ne.0) then
   call warning(whoami, 'Error encountered while reading file "'//fname(1:flen)//'.',-1)
  else
!********* read scalar *******
! note: will implement error control, in case file is corrupted
  if (nx.ne.mx.or.ny.ne.my.or.nz.ne.mz) then
!write(0,*) nx,mx,ny,my,nz,mz ! aa
    call warning(whoami, 'Mismatch in the domain size. Abort', -1)
  else
   if (q3d) then
    if (qbin) then ; read(fid, IOSTAT=ioerr) q ; else ; read(fid,*, IOSTAT=ioerr) q ; endif
   else
    if (qbin) then ; read(fid, IOSTAT=ioerr) q(:,:,2) ; else ; read(fid,*, IOSTAT=ioerr) q(:,:,2) ; endif
   endif
  endif ! nx==mx
!
  if (ioerr.ne.0) call warning(whoami, 'Error encountered while reading file "'//fname(1:flen)//'.',-1)
!
  endif ! ioerr
!
  call files_close(fid)
!
  endif ! me.eq.0
!
  if (.not.fatal_warning()) then
   if (bin) then ; call message(whoami, 'Binary scalar file "'//fname(1:flen)//'" read.')
   else ; call message(whoami, 'ASCII scalar file "'//fname(1:flen)//'" read.')
   endif
!
  endif ! fatal_warning
!
  end subroutine chest_read_scalar
