
interface 
 function erfsun(a)
 real*8 :: a, erfsun
 end function erfsun
end interface

! real*8 :: erfsun
!
 real*8 :: a=0.12345d0, b, c
!
! do j=1,10
! do i=1,10000000
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
