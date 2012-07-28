
  implicit none
  float :: x1(66), y1(66), z1(66)
  float :: x2(66), y2(66), z2(66)
  float :: xbc0(66,66), xbc1(66,66)
  float :: ybc0(66,66), ybc1(66,66)
  float :: zbc0(66,66), zbc1(66,66)
  float :: rhs(66,66,66), u(66,66,66)
  float, parameter :: err = 1d-12
  float :: a, b, l, z, pi, delta
  float :: fx, fy, fz
  integer :: i, j, k, n
!
!
! xgrid:
  a=0d0; b=1d0; z=1.05d0; n=66-1
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
  x2(2:66)=0.5d0*(x1(1:66-1)+x1(2:66)); x2(1)=x1(1)-0.5d0*(x1(2)-x1(1)); 
!
! ygrid:
  a=0d0; b=2d0; z=1.05d0; n=66-1
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
  y2(2:66)=0.5d0*(y1(1:66-1)+y1(2:66)); y2(1)=y1(1)-0.5d0*(y1(2)-y1(1)); 
!
! zgrid:
  a=0d0; b=1.5d0; z=1.05d0; n=66-1
  l=b-a 
  if (ABS(z-1.0d0).lt.ERR) then 
   delta=l/(n-1) 
  else 
   delta=l*(1.0d0-z)/(1.0d0-(z**(n-1)))
  endif
!
  z1(1)=a
  do i=2,n+1
   z1(i)=z1(i-1)+delta;
   delta=delta*z;
  enddo
!
  open(1,file='zg_test.dat',form='formatted', status='unknown')
  write(1,'(F30.15)') z1 
  close(1)
!
  z2(2:66)=0.5d0*(z1(1:66-1)+z1(2:66)); z2(1)=z1(1)-0.5d0*(z1(2)-z1(1)); 
!
! data arrays
!
 pi=atan(1d0)*4d0
 do i=1,66 ; do j=1, 66; do k=1, 66
  fx=1d0-exp(-4d0*x2(i))
  fy=cos(2d0*pi*1*y2(j))
  fz=sin(2d0*pi*2*z2(k))
  u(i,j,k) = fx*fy*fz
  rhs(i,j,k) = 4d0**2 * (fx-1d0) * fy * fz - 4d0*pi*pi*(1*1+2*2)*u(i,j,k)
 enddo;       enddo;       enddo
!
! write datafiles in chest format
!
 open(1,file='uexact_test.dat',form='formatted', status='unknown')
 write(1,'(3I5)') 66,66,66
 write(1,'(66F30.15)') u 
! write(1) 66,66,66
! write(1) u 
 close(1)
 open(1,file='rhs_test.dat',form='formatted', status='unknown')
 write(1,'(3I5)') 66,66,66
 write(1,'(66F30.15)') rhs 
! write(1) 66,66,66
! write(1) rhs 
 close(1)
! bc arrays
 a=(1d0-exp (-4d0*0))
 b=(1d0-exp (-4d0*1))
 do k=1,66 ; fz=sin(2d0*pi*2*z2(k))
  do j=1,66 ; fy=cos(2d0*pi*1*y2(j))
   xbc0(j,k) = a * fy * fz
   xbc1(j,k) = b * fy * fz
  enddo
 enddo
!
 a=cos(2d0*pi*1*0)
 b=cos(2d0*pi*1*2)
 do k=1,66 ; fz=sin(2d0*pi*2*z2(k))
  do i=1,66 ; fx=(1d0-exp (-4d0*x2(i)))
   ybc0(i,k) = fx * a * fz
   ybc1(i,k) = fx * b * fz
  enddo
 enddo
!
 a=sin(2d0*pi*2*0)
 b=sin(2d0*pi*2*1.5)
 do j=1,66 ; fy=cos(2d0*pi*1*y2(j))
  do i=1,66 ; fx=(1d0-exp (-4d0*x2(i)))
   zbc0(i,j) = fx * fy * a
   zbc1(i,j) = fx * fy * b
  enddo
 enddo
!
 open(10,file='xbc0_test.dat',form='formatted', status='unknown')
 open(11,file='xbc1_test.dat',form='formatted', status='unknown')
 open(12,file='ybc0_test.dat',form='formatted', status='unknown')
 open(13,file='ybc1_test.dat',form='formatted', status='unknown')
 open(14,file='zbc0_test.dat',form='formatted', status='unknown')
 open(15,file='zbc1_test.dat',form='formatted', status='unknown')
!
 write(10,'(F30.15)') xbc0(2:66-1,2:66-1); close(10)
 write(11,'(F30.15)') xbc1(2:66-1,2:66-1); close(11)
 write(12,'(F30.15)') ybc0(2:66-1,2:66-1); close(12)
 write(13,'(F30.15)') ybc1(2:66-1,2:66-1); close(13)
 write(14,'(F30.15)') zbc0(2:66-1,2:66-1); close(14)
 write(15,'(F30.15)') zbc1(2:66-1,2:66-1); close(15)
!
end 
