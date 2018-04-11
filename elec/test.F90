
program test

 use vars

 int :: n=100
 int i
 float :: errx, erry, errz

 __FARR(mem,:,:)
 __FARR(f,:);
 __FARR(r,:);
 __FPAR(h,0.00001d0)
 __FPAR(oo2h,1d0/(2d0*h))
!
!
 __FARR(elsr_dfd,:,:)
 __FARR(elsr_dx_fd,:)
 __FARR(elsr_dy_fd,:)
 __FARR(elsr_dz_fd,:)
 
 __OUT(' Running tests')
 __OUT(' Calling setup')
 call setup()
 
 __OUT(' Check fshortp function: ')
 __ALLOC(mem(n,2))
 r=>mem(:,1);
 f=>mem(:,2);
 do i=1,n
  r(i)=1d0*i/n;
  f(i)=fshortp(r(i))
 enddo
 
 call write_field(mem,'fshortp.txt',ascii,2*n)
 
 
 __OUT(' FD test of short-range potential :' )
 __ALLOC(elsr_dfd(npt,3))
 elsr_dfd=0d0
!
 elsr_dx_fd=>elsr_dfd(:,1)
 elsr_dy_fd=>elsr_dfd(:,2)
 elsr_dz_fd=>elsr_dfd(:,3)
!
 __OUT(' Computing analytical gradients');
 call filt3p()
 __OUT(' Computing finite-difference gradients');
 do i=1,npt
! x-grad :
  x(i)=x(i)+h ; 
  call filt3p()
  elsr_dx_fd(i)=elsr ; 

  x(i)=x(i)-h-h ;
  call filt3p()
  __INCR(elsr_dx_fd(i),-elsr)
 __SCALE(elsr_dx_fd(i),oo2h)

! y-grad : 
  y(i)=y(i)+h ; 
  call filt3p()
  elsr_dy_fd(i)=elsr ; 

  y(i)=y(i)-h-h ;
  call filt3p()
  __INCR(elsr_dy_fd(i),-elsr)
 __SCALE(elsr_dy_fd(i),oo2h)
! z-grad : 
  z(i)=z(i)+h ; 
  call filt3p()
  elsr_dz_fd(i)=elsr ; 

  z(i)=z(i)-h-h ;
  call filt3p()
  __INCR(elsr_dz_fd(i),-elsr)
 __SCALE(elsr_dz_fd(i),oo2h)
 enddo ! npt

 errx=sum(abs(elsr_dx-elsr_dx_fd));
 erry=sum(abs(elsr_dy-elsr_dy_fd));
 errz=sum(abs(elsr_dz-elsr_dz_fd));

 __OUT(' Total absolute x-grad error : ', errx )
 __OUT(' Total absolute y-grad error : ', erry )
 __OUT(' Total absolute z-grad error : ', errz )
 if ((errx+erry+errz)<h) then
  __OUT( ' Total absolute grad error < h (=',h,'). PASSED.')
 else
  __OUT( ' Total absolute grad error >= h (=',h,'). FAILED.')
 endif

end
