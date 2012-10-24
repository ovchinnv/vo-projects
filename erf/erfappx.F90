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
! 12.20.11: did not find using 5th order approximation faster than the 7th order one
! this may be due to the cost of the exponential function (as well as other operations within the loop)
! may be a good idea to approximate the exponential function
!DEC$ ATTRIBUTES FORCEINLINE :: erfo7
!DEC$ ATTRIBUTES FORCEINLINE :: erfo5
! approximations from Abramowits & Stegun (also see Wiki)
function erfo7(y)
implicit none
real*8, parameter :: a1=0.0705230784d0,a2=0.0422820123d0,a3=0.0092705272d0,a4=0.0001520143d0,a5=0.0002765672d0,a6=0.0000430638d0,&
& one=1d0, zero=0d0
real*8 :: erfo7, y, x, x2, x3
integer :: isgn
isgn=sign(one,y); x=isgn*y;
x2=x*x; x3=x2*x
erfo7 = one + a1 * x + a2 * x2 + a3 * x3 + a4 * x2 * x2 + a5 * x2 * x3 + a6 * x3 * x3;
erfo7=one/erfo7
erfo7=erfo7*erfo7; erfo7=erfo7*erfo7; erfo7=erfo7*erfo7; erfo7=erfo7*erfo7; ! 16th power
erfo7=(one-erfo7)*isgn
end function erfo7
!
function erfo5(y)
implicit none
real*8, parameter :: p=0.47047d0, a1=0.3480242d0, a2=-0.0958798d0, a3=0.7478556d0, one=1d0, zero=0d0
real*8 :: erfo5, y, x, t, t2
integer :: isgn
isgn=sign(one,y); x=isgn*y;
t=one/(one+p*x); t2=t*t
erfo5 = (one - (a1*t + a2*t2 + a3*t*t2)*exp(-x*x))*isgn
end function erfo5
!
