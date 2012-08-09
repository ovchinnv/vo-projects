/*#define __WRN(__WHO,__MSG) write(0,*) 'WARNING FROM: ',__WHO,': ',__MSG*/
/*#define __PRINT(__MSG) write(0,'(A)') __MSG*/
/*#define __PRINT(__MSG) call plainmessage(__MSG)*/
/*#define __PRINTL(__MSG,__LEVEL) call plainmessage(__MSG,__LEVEL)*/
/*COORDINATES AND MASSES:*/
/*#define __INDX(__STR, __STRLEN, __TEST, __TESTLEN)  index(__STR(1:min(__STRLEN,len(__STR))),__TEST(1:min(__TESTLEN,len(__TEST))))*/



/*
#ifdef __PRINT
#undef __PRINT
#endif
#define __PRINT(__WHAT) call plainmessage(__WHAT)
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
! subroutines to compute various statistical quantities
module stats
 use constants
 contains
 function calcKE(v, m)
 implicit none
 integer :: natom, ndim
 real*8 :: v(:,:), m(:)
 real*8, dimension (2) :: calcKE
! assume that the dimensions are correct for multiplication to make sense
 natom=size(v,2)
 ndim=size(v,1)
!
 calcKE(1)=0.5d0*sum(matmul(v**2,m))
 calcKE(2)=calcKE(1)*2/ndim/kboltzmann/natom
!
 end function calcKE
end module stats
