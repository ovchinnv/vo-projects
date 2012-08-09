/*#define __WRN(__WHO,__MSG) write(0,*) 'WARNING FROM: ',__WHO,': ',__MSG*/
/*#define __PRINT(__MSG) write(0,'(A)') __MSG*/
/*#define __PRINT(__MSG) call plainmessage(__MSG)*/
/*#define __PRINTL(__MSG,__LEVEL) call plainmessage(__MSG,__LEVEL)*/
/*COORDINATES AND MASSES:*/
/*#define __INDX(__STR, __STRLEN, __TEST, __TESTLEN)  index(__STR(1:min(__STRLEN,len(__STR))),__TEST(1:min(__TESTLEN,len(__TEST))))*/
/*#define __WRN(__WHO,__MSG) write(0,*) 'WARNING FROM: ',__WHO,': ',__MSG*/
/*#define __PRINT(__MSG) write(0,'(A)') __MSG*/
/*#define __PRINT(__MSG) call plainmessage(__MSG)*/
/*#define __PRINTL(__MSG,__LEVEL) call plainmessage(__MSG,__LEVEL)*/
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
   subroutine writeblock(f,q,nx,ny,nz,qb,q3)
   implicit none
   integer :: f,i,j,k,nx,ny,nz
   real*8 :: q(nx,ny,nz)
   logical :: qb, q3
!
   if (qb) then
    write(f) 0,0,0,0 ! plot 3d format specifies four reals ( Mach number, angle of attack, Reynolds number, time )
    if (q3) then
       write(f) (((real(q(i,j,k)),i=1,nx),j=1,ny),k=1,nz)
    else
     write(f) ((( real(q(i,j,k)),i=1,nx) ,j=1,ny), k=2,2)
    endif
   else
    write(f,*) 0,0,0,0 ! see above
    if (q3) then
     write(f,*) (((q(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    else
     write(f,*) (((q(i,j,k),i=1,nx),j=1,ny),k=2,2) ! write the single inner point
    endif
   endif
   end subroutine writeblock
!
!*******************************************************************************************************************
   subroutine readblock(f,q,nx,ny,nz,qb,q3)
   implicit none
   integer :: f,i,j,k,nx,ny,nz
   real*8 :: q(nx,ny,nz)
   real*4 :: qsngl(nx,ny,nz)
   real*8 :: mach, alpha, Re, t
   logical :: qb, q3
!
   if (qb) then
    read(f) mach, alpha, Re, t ! plot 3d format specifies four reals ( Mach number, angle of attack, Reynolds number, time )
    if (q3) then
       read(f) (((qsngl(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    else
     read(f) (((qsngl(i,j,k),i=1,nx) ,j=1,ny), k=2,2)
    endif
    q=qsngl
   else
    read(f,*) mach, alpha, Re, t ! see above
    if (q3) then
     read(f,*) (((q(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    else
     read(f,*) (((q(i,j,k),i=1,nx),j=1,ny),k=2,2) ! read the single inner point
    endif
   endif
   end subroutine readblock
!
!*******************************************************************************************************************
! Subroutine to write five (four in 2D) quantities at once; this is the Plot3D standard
   subroutine writeblock5(f,q1,q2,q3,q4,q5,nx,ny,nz,qb,q3d)
   implicit none
   integer :: f,i,j,k,nx,ny,nz
   real*8, dimension(nx,ny,nz) :: q1,q2,q3,q4,q5
   logical :: qb, q3d
!
   if (qb) then
    write(f) 0,0,0,0 ! plot 3d format specifies four reals ( Mach number, angle of attack, Reynolds number, time )
    if (q3d) then
       write(f) &
& (((real(q1(i,j,k)),i=1,nx),j=1,ny),k=1,nz),&
& (((real(q2(i,j,k)),i=1,nx),j=1,ny),k=1,nz),&
& (((real(q3(i,j,k)),i=1,nx),j=1,ny),k=1,nz),&
& (((real(q4(i,j,k)),i=1,nx),j=1,ny),k=1,nz),&
& (((real(q5(i,j,k)),i=1,nx),j=1,ny),k=1,nz)
    else
     write(f) &
& ((( real(q1(i,j,k)),i=1,nx) ,j=1,ny),k=2,2),&
& ((( real(q2(i,j,k)),i=1,nx) ,j=1,ny),k=2,2),&
& ((( real(q3(i,j,k)),i=1,nx) ,j=1,ny),k=2,2),&
& ((( real(q4(i,j,k)),i=1,nx) ,j=1,ny),k=2,2)
    endif
   else
    write(f,*) 0,0,0,0 ! see above
    if (q3d) then
     write(f,*) &
& (((real(q1(i,j,k)),i=1,nx),j=1,ny),k=1,nz),&
& (((real(q2(i,j,k)),i=1,nx),j=1,ny),k=1,nz),&
& (((real(q3(i,j,k)),i=1,nx),j=1,ny),k=1,nz),&
& (((real(q4(i,j,k)),i=1,nx),j=1,ny),k=1,nz),&
& (((real(q5(i,j,k)),i=1,nx),j=1,ny),k=1,nz)
    else
     write(f,*) &
& (((real(q1(i,j,k)),i=1,nx),j=1,ny),k=2,2),&
& (((real(q2(i,j,k)),i=1,nx),j=1,ny),k=2,2),&
& (((real(q3(i,j,k)),i=1,nx),j=1,ny),k=2,2),&
& (((real(q4(i,j,k)),i=1,nx),j=1,ny),k=2,2)
    endif
   endif
   end subroutine writeblock5
!
