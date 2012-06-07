
  implicit none
  float :: x1(258), y1(258)
  float :: x2(258), y2(258)
  float :: xbc0(258), xbc1(258)
  float :: ybc0(258), ybc1(258)
  float :: rhs(258,258), u(258,258)
  float, parameter :: err = 1d-12
  float :: a, b, l, z, pi, delta
  integer :: i, j, n
!
! xgrid:
  a=0d0; b=1d0; z=1.0d0; n=258-1
  l=b-a
  if (ABS(z-1.0d0).lt.ERR) then 
   delta=l/(n-1) 
  else 
   delta=l*(1.0d0-z)/(1.0d0-(z**(n-1)))
  endif
!
  x1(1)=a
  do i=2,n+1
   x1(i)=x1(i-1)+delta;
   delta=delta*z;
  enddo
!
  open(1,file='xg_test.dat',form='formatted', status='unknown')
  write(1,'(F30.15)') x1 
  close(1)
!
  x2(2:258)=0.5d0*(x1(1:258-1)+x1(2:258)); x2(1)=x1(1)-0.5d0*(x1(2)-x1(1)); 
!
! ygrid:
  a=0d0; b=2d0; z=1.0d0; n=258-1
  l=b-a 
  if (ABS(z-1.0d0).lt.ERR) then 
   delta=l/(n-1) 
  else 
   delta=l*(1.0d0-z)/(1.0d0-(z**(n-1)))
  endif
!
  y1(1)=a
  do i=2,n+1
   y1(i)=y1(i-1)+delta;
   delta=delta*z;
  enddo
!
  open(1,file='yg_test.dat',form='formatted', status='unknown')
  write(1,'(F30.15)') y1 
  close(1)
!
  y2(2:258)=0.5d0*(y1(1:258-1)+y1(2:258)); y2(1)=y1(1)-0.5d0*(y1(2)-y1(1)); 
!
! data arrays
!
 pi=atan(1d0)*4d0
 do i=1,258 ; do j=1, 258
  u(i,j) = (1-exp (-5d0*x2(i)))*cos(2*pi*1*y2(j))
  rhs(i,j) = -5d0**2 * exp(-5d0*x2(i))*cos(2*pi*1*y2(j)) - (2*pi*1)**2 * u(i,j)
 enddo; enddo
!
! write datafiles in chest format
!
 open(1,file='uexact_test.dat',form='formatted', status='unknown')
 write(1,'(3I5)') 258,258
 write(1,'(258F30.15)') u 
! write(1) 258,258
! write(1) u 
 close(1)
 open(1,file='rhs_test.dat',form='formatted', status='unknown')
 write(1,'(3I5)') 258,258
 write(1,'(258F30.15)') rhs 
! write(1) 258,258
! write(1) rhs 
 close(1)
! bc arrays
 xbc0= (1-exp (-5d0*0))  *cos(2*pi*1*y2(:))
 xbc1= (1-exp (-5d0*1))  *cos(2*pi*1*y2(:))
 ybc0= (1-exp (-5d0*x2(:)))  *cos(2*pi*1*0)
 ybc1= (1-exp (-5d0*x2(:)))  *cos(2*pi*1*2)
!
 open(10,file='xbc0_test.dat',form='formatted', status='unknown')
 open(11,file='xbc1_test.dat',form='formatted', status='unknown')
 open(12,file='ybc0_test.dat',form='formatted', status='unknown')
 open(13,file='ybc1_test.dat',form='formatted', status='unknown')
!
 write(10,'(F30.15)') xbc0(2:258-1); close(10)
 write(11,'(F30.15)') xbc1(2:258-1); close(11)
 write(12,'(F30.15)') ybc0(2:258-1); close(12)
 write(13,'(F30.15)') ybc1(2:258-1); close(13)
!
end 
