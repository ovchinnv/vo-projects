

interface erf_sun
 function erfsun_r(a)
 real*8 :: a, erfsun_r
 end function erfsun_r
end interface


!use erf_ext

  real*8 :: a=0.12345d0, b, c, d, e, f

  b=erf_sun(a)
!  b=erfsun(a)
!  c=erf_c(a)
  d=erf(a)

!  e=erfo7(a)  ;  f=erfo5(a)  
  
!
!  write(0,*) b,c,e,f,d
  write(0,*) b,d
!
! 
 end
