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
!module plot3Dio
! use output
!
! contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine plot3Dwrite_grid(filename,x,y,z,nx,ny,nz,nblock,qbin)
  use gridsize, only : me, communicator
  use output
  use files
  use parser
  implicit none
!
  integer :: fid=-1
  real*8 :: x(*),y(*),z(*)
  integer :: nx(*), ny(*), nz(*) ! have dimensions of nblock[if present] or 1, whichever is greater
! there may be problems with optional arguments here
! integer, optional :: nblock
! logical, optional :: qbin
  integer :: nblock
  logical :: qbin
  character(len=11) :: fmt
!
  integer :: nb
  logical :: bin
  logical :: q3d
!
  integer :: i, j, k, m, ioff, joff, koff
!
  character(len=*) :: filename
  character(len=len(filename)) :: fname
  character(len=16) , parameter :: whoami='PLOT3DWRITE_GRID'
  integer :: flen
!
  if (me.gt.0) return ! may need better parallel treatment
!
! if ( present(nblock) ) then ; nb=nblock ; else ; nb=1 ; endif;
! if ( present(qbin ) ) then ; bin=qbin ; else ; bin=.false. ; endif;
  nb=nblock
  bin=qbin
  q3d=(nz(1).gt.3)
!
  if (bin) then; fmt='UNFORMATTED'; else ; fmt='FORMATTED'; endif
!
  fname=filename
  call adjustleft(fname)
  flen=len_trim(fname)
  if (flen.gt.0) then
   if (me.eq.0) then
    call files_open(fid, name_=fname(1:flen), form_=fmt, action_='WRITE')
    if (fid.lt.0) then
     call warning(whoami, 'Cannot open grid file. Abort.',-1)
    endif
   endif ! me
   if (fatal_warning()) return
!
  else
   call warning(whoami, 'Grid file name not specified. Abort.',0)
   return
  endif
!
  if (bin) then
   if (nb.gt.1) write(fid) nb
   if (q3d) then;write(fid) ( nx(m), ny(m), nz(m), m = 1, nb ) ;
            else;write(fid) ( nx(m), ny(m), m = 1, nb ) ; endif
   ioff=0; joff=0; koff=0;
   do m = 1, nblock
    if (q3d) then
     write(fid) (((real(x(i+ioff)), i=1,nx(m)), j=1,ny(m)), k=1, nz(m)), &
& (((real(y(j+joff)), i=1,nx(m)), j=1,ny(m)), k=1, nz(m)), &
& (((real(z(k+koff)), i=1,nx(m)), j=1,ny(m)), k=1, nz(m))
     koff=koff+nz(m);
    else
     write(fid) ((real(x(i+ioff)), i=1,nx(m)), j=1,ny(m)), &
& ((real(y(j+joff)), i=1,nx(m)), j=1,ny(m))
    endif
!
    ioff=ioff+nx(m)
    joff=joff+ny(m)
!
   enddo ! blocks
!
  else ! bin
!aa
   if (nb.gt.1) write(fid,*) nb
   if (q3d) then;write(fid,*) ( nx(m), ny(m), nz(m), m = 1, nb ) ;
            else;write(fid,*) ( nx(m), ny(m), m = 1, nb ) ; endif
   ioff=0; joff=0; koff=0;
   do m = 1, nb
    if (q3d) then
     write(fid,*) (((x(i+ioff), i=1,nx(m)), j=1,ny(m)), k=1, nz(m)), &
& (((y(j+joff), i=1,nx(m)), j=1,ny(m)), k=1, nz(m)), &
& (((z(k+koff), i=1,nx(m)), j=1,ny(m)), k=1, nz(m))
     koff=koff+nz(m);
    else
     write(fid,*) ((x(i+ioff), i=1,nx(m)), j=1,ny(m)), &
& ((y(j+joff), i=1,nx(m)), j=1,ny(m))
    endif
!
    ioff=ioff+nx(m)
    joff=joff+ny(m)
!
   enddo
  endif ! qbin
!
  call files_close(fid)
!
  if (bin) then ; call message(whoami, 'Binary grid file "'//fname(1:flen)//'" written.')
  else ; call message(whoami, 'ASCII grid file "'//fname(1:flen)//'" written.')
  endif
!
  end subroutine plot3Dwrite_grid
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine plot3Dwrite_scalar(filename,q,nx,ny,nz,nblock,qbin) ! scalar function
  use output
  use files
  use gridsize, only : me, communicator
  use parser
  implicit none
!
  real*8 :: q(*)
  integer :: nx(*), ny(*), nz(*) ! have dimensions of nblock[if present] or 1, whichever is greater
! integer, optional :: nblock
! logical, optional :: qbin
  integer :: nblock
  logical :: qbin
  character(len=11) :: fmt
  character(len=18) , parameter :: whoami='PLOT3DWRITE_SCALAR'
!
  integer :: nb, offset
  logical :: bin
  logical :: q3d
!
  integer :: i, j, k, m, ioff, joff, koff
!
  character(len=*) :: filename
  character(len=len(filename)) :: fname
  integer :: flen, fid=-1
!
! if ( present(qbin ) ) then ; bin=qbin ; else ; bin=.true. ; endif;
  bin=qbin
  if (bin) then; fmt='UNFORMATTED'; else ; fmt='FORMATTED'; endif
!
  fname=filename
  call adjustleft(fname)
  flen=len_trim(fname)
  if (flen.gt.0) then
   if (me.eq.0) then
    call files_open(fid, name_=fname(1:flen), form_=fmt, action_='WRITE')
    if (fid.lt.0) then
     call warning(whoami, 'Cannot open output file. Abort.',-1)
    endif
   endif ! me
   if (fatal_warning()) return
  else
   call warning(whoami, 'Output file name not specified. Abort.',0)
   return
  endif
!
  if (me.gt.0) return ! may change parallel behavior in the future
!
! if ( present(nblock) ) then ; nb=nblock ; else ; nb=1 ; endif;
  nb=nblock
  if (nz(1).gt.3) then ; q3d=.true. ; else ; q3d=.false.; endif;
!************** header *******
!
  if (bin) then
   if (nb.gt.1) write(fid) nb
   if (q3d) then
    write(fid) ( nx(m), ny(m), nz(m), m = 1, nb )
   else
    write(fid) ( nx(m), ny(m), 1, m = 1, nb )
   endif ! q3d
  else ! qbin
   if (nb.gt.1) write(fid,*) nb
   if (q3d) then
    write(fid,*) ( nx(m), ny(m), nz(m), m = 1, nb )
   else
    write(fid,*) ( nx(m), ny(m), m = 1, nb )
   endif ! q3d
  endif
!
!********* write scalar *******
  offset=0
  do m=1, nblock
   call writeblock(fid,q(offset+1),nx(m),ny(m),nz(m),bin,q3d) ! trick fortran into redimensioning q
   offset=offset+nx(m)*ny(m)*nz(m)
  enddo
!
  call files_close(fid)
!
  if (bin) then ; call message(whoami, 'Binary scalar file "'//fname(1:flen)//'" written.')
  else ; call message(whoami, 'ASCII scalar file "'//fname(1:flen)//'" written.')
  endif
!
  end subroutine plot3Dwrite_scalar
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine plot3Dread_scalar(filename,q,nx,ny,nz,nblock,qbin) ! scalar function
  use output
  use files
  use gridsize, only : me, communicator, ncpu
  use parser
  implicit none
!
  real*8 :: q(*)
  integer :: nx(*), ny(*), nz(*) ! have dimensions of nblock[if present] or 1, whichever is greater
  integer, allocatable, dimension(:) :: mx, my, mz
! integer, optional :: nblock
! logical, optional :: qbin
  integer :: nblock
  logical :: qbin
  character(len=11) :: fmt
  character(len=17) , parameter :: whoami='PLOT3DREAD_SCALAR'
!
  integer :: nb, offset, bug
  logical :: bin
  logical :: q3d
!
  integer :: i, j, k, m, ioff, joff, koff
!
  character(len=*) :: filename
  character(len=len(filename)) :: fname
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
  if (me.eq.0) then
!***************** ONLY ROOT READS, THEN BROADCASTS TO OTHER NODES, IF NEEDED ****
  if (nz(1).gt.3) then ; q3d=.true. ; else ; q3d=.false.; endif;
!************** header *******
!
  if (bin) then
   if (nblock.gt.1) read(fid) nb
   if (nb.ne.nblock) &
& call warning(whoami, 'Mismatch in the number of grid blocks,',-1)
   allocate(mx(nblock), my(nblock), mz(nblock))
!
   if (q3d) then
    read(fid) ( mx(m), my(m), mz(m), m = 1, nb )
   else
    read(fid) ( mx(m), my(m), m = 1, nb )
   endif ! q3d
  else ! qbin
   if (nb.gt.1) read(fid,*) nb
   if (nb.ne.nblock) &
& call warning(whoami, 'Mismatch in the number of grid blocks.',-1)
   allocate(mx(nblock), my(nblock), mz(nblock))
   if (q3d) then
    read(fid,*) ( mx(m), my(m), mz(m), m = 1, nb )
   else
    read(fid,*) ( mx(m), my(m), m = 1, nb )
   endif ! q3d
  endif
!
!********* read scalar *******
  offset=0
  do m=1, nblock
   if (nx(m).ne.mx(m).or.ny(m).ne.my(m).or.nz(m).ne.mz(m)) then
    call warning(whoami, 'Mismatch in the domain size. Abort', -1)
    exit
   endif
   call readblock(fid,q(offset+1),nx(m),ny(m),nz(m),bin,q3d) ! trick fortran into redimensioning q
   offset=offset+nx(m)*ny(m)*nz(m)
  enddo
!
  call files_close(fid)
!
  if (allocated(mx)) deallocate(mx)
  if (allocated(my)) deallocate(my)
  if (allocated(mz)) deallocate(mz)
!
  endif ! me.eq.0 (Root reads)
!
  if (.not.fatal_warning()) then
   if (bin) then ; call message(whoami, 'Binary scalar file "'//fname(1:flen)//'" read.')
   else ; call message(whoami, 'ASCII scalar file "'//fname(1:flen)//'" read.')
   endif
!
!
  endif ! fatal_warning
!
  end subroutine plot3Dread_scalar
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine plot3Dwrite_solution(filename,q1,q2,q3,q4,q5,nx,ny,nz,nblock,qbin)
! write 4 or five scalar function at once, which represent the solution in the standard plot3D format
  use output
  use files
  use gridsize, only : me, communicator
  use parser
  implicit none
!
  real*8, dimension(*) :: q1, q2, q3, q4, q5
  integer, dimension(*) :: nx, ny, nz ! have dimensions of nblock[if present] or 1, whichever is greater
! integer, optional :: nblock
! logical, optional :: qbin
  integer :: nblock
  logical :: qbin
  character(len=11) :: fmt
  character(len=20) , parameter :: whoami='PLOT3DWRITE_SOLUTION'
!
  integer :: nb, offset
  logical :: bin
  logical :: q3d
!
  integer :: i, j, k, m, ioff, joff, koff
!
  character(len=*) :: filename
  character(len=len(filename)) :: fname
  integer :: flen, fid=-1
!
! if ( present(qbin ) ) then ; bin=qbin ; else ; bin=.true. ; endif;
  bin=qbin
  if (bin) then; fmt='UNFORMATTED'; else ; fmt='FORMATTED'; endif
!
  fname=filename
  call adjustleft(fname)
  flen=len_trim(fname)
  if (flen.gt.0) then
   if (me.eq.0) then
    call files_open(fid, name_=fname(1:flen), form_=fmt, action_='WRITE')
    if (fid.lt.0) then
     call warning(whoami, 'Cannot open output file. Abort.',-1)
    endif
   endif ! me
   if (fatal_warning()) return
  else
   call warning(whoami, 'Output file name not specified. Abort.',0)
   return
  endif
!
  if (me.gt.0) return ! may change parallel behavior in the future
!
! if ( present(nblock) ) then ; nb=nblock ; else ; nb=1 ; endif;
  nb=nblock
  if (nz(1).gt.3) then ; q3d=.true. ; else ; q3d=.false.; endif;
!************** header *******
!
  if (bin) then
   if (nb.gt.1) write(fid) nb
   if (q3d) then
    write(fid) ( nx(m), ny(m), nz(m), m = 1, nb )
   else
    write(fid) ( nx(m), ny(m), m = 1, nb )
   endif ! q3d
  else ! qbin
   if (nb.gt.1) write(fid,*) nb
   if (q3d) then
    write(fid,*) ( nx(m), ny(m), nz(m), m = 1, nb )
   else
    write(fid,*) ( nx(m), ny(m), m = 1, nb )
   endif ! q3d
  endif
!
!********* write scalar *******
  offset=0
  do m=1, nblock
   call writeblock5(fid,q1(offset+1),q2(offset+1),q3(offset+1),q4(offset+1),q5(offset+1),nx(m),ny(m),nz(m),bin,q3d) ! trick fortran into redimensioning q
   offset=offset+nx(m)*ny(m)*nz(m)
  enddo
!
  call files_close(fid)
!
  if (bin) then ; call message(whoami, 'Binary solution file "'//fname(1:flen)//'" written.')
  else ; call message(whoami, 'ASCII solution file "'//fname(1:flen)//'" written.')
  endif
!
  end subroutine plot3Dwrite_solution
!
!end module plot3Dio
