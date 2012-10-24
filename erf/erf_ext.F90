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
! NOTE: I have not been able to interface C and Fortran without the "value" parameter in the interface, which requires 2003 standard
module erf_ext
!
use, intrinsic :: iso_c_binding
implicit none
!
private
public erf_sun
public erf_c
!
interface
 function derf_sun(a) bind(c,NAME='erfsun_v_')
 import
 real(c_double), intent(in), value :: a
 real(c_double) :: derf_sun
 end function derf_sun
!
! function from C stdlib
!
 function derf_c(x) bind(c,NAME='erf')
 import
 real(c_double), intent(in), value :: x
 real(c_double) :: derf_c
 end function derf_c
end interface
interface erf_c
 module procedure erf_double
 module procedure erf_single
end interface
!
!
interface erf_sun
 module procedure erfsun_double
 module procedure erfsun_single
end interface
!
contains
!_________________________WRAPPER FUNCTIONS ___________________
!
function erf_double(x)
real*8, intent(in) :: x
real*8 :: erf_double
real(c_double) :: y
y=x ; erf_double=derf_c(y)
end function erf_double
!
function erf_single(x)
real*4, intent(in) :: x
real*4 :: erf_single
real(c_double) :: y
y=x ; erf_single=derf_c(y)
end function erf_single
!
function erfsun_double(x)
real*8, intent(in) :: x
real*8 :: erfsun_double
real(c_double) :: y
y=x ; erfsun_double=derf_sun(y)
end function erfsun_double
!
function erfsun_single(x)
real*4, intent(in) :: x
real*4 :: erfsun_single
real(c_double) :: y
y=x ; erfsun_single=derf_sun(y)
end function erfsun_single
!
!
end module erf_ext
