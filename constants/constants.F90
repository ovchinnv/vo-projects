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
module constants
!
 public
!
 real*8, parameter :: pi = 3.141592653589793d0
 real*8, parameter :: twopi = 3.141592653589793d0*2d0
!
 real*8, parameter :: kboltzmann = 1.987191d-3 ! Boltzmann constant
 real*8, parameter :: fs_per_akma = 48.8882129d0 ! factor to convert between AKMA time units and picoseconds
 real*8, parameter :: akma_per_fs=1d0/fs_per_akma
 real*8, parameter :: zero=0d0
 real*8, parameter :: one=1d0
 real*8, parameter :: two=2d0
 real*8, parameter :: three=3d0
 real*8, parameter :: nine=9d0
 real*8, parameter :: twentyseven=27d0
 real*8, parameter :: half=0.5d0
 real*8, parameter :: quarter=0.25d0
 real*8, parameter :: third=1d0/3d0
 real*8, parameter :: unknownf=1d0*(not(ishftc(1,-1))) ! smallest representable integer
!
! real*8, parameter :: DTOL = 1.0e-14
! real*8, parameter :: FTOL = 1.0e-7
 real*8, save, private :: TOL = -1 ! error tolerance
!**********************************************************
 contains
  function ERRTOL()
   implicit none
   real*8 :: a, b, ERRTOL
   if (TOL.lt.0) then ! compute machine precision
    a=1.0; b=1.0;
    do
     if ( a - a / b .ne. a ) then ; b = b * 10 ; else ; TOL = 50 / b ; exit ; endif
    enddo
   endif
   ERRTOL=TOL
  end function ERRTOL
!
end module constants
