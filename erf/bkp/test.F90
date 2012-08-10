
use, intrinsic :: iso_c_binding

interface 
 function erfsun(a)
 import
 real(c_double) :: a, erfsun
 end function erfsun
end interface

! real*8 :: erfsun
!
 real*8 :: a=0.12345d0, b, c
 real( c_double ) :: a_
!
! do j=1,10
! do i=1,10000000
   a_=a
  b=erfsun(a)
!  a=100.
  c=erf(a)
!  c=1./a
!  c=sqrt(a)
!  c=exp(a*a)
!  a=a+0.00001
  write(0,*) b,c
!  a=a+0.0000001d0
! enddo 
! enddo
!
! write(0,*) a, b !, c
! 
 end
