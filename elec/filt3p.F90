
subroutine filt3p()
 use vars

 int :: nimg=1; ! number of images to consider
 float :: oos=1/spt;
 float :: s2=spt*spt;
 float r2q(npt);
 float dxq(npt);
 float dyq(npt);
 float dzq(npt);

 float xmin, xmax, ymin, ymax, zmin, zmax
 int imin, imax, jmin, jmax, kmin, kmax
 int ii, jj, kk, img, jmg, kmg
 int i, j
 float r, r2
 float dx2, dy2, dz2, odx, ody, odz
 float dx, dy, dz, pre

 odx=1d0/(xx(2)-xx(1)) ; !assume uniform
 ody=1d0/(yy(2)-yy(1))
 odz=1d0/(zz(2)-zz(1))
!
 rho=0d0 ! initialize
 elsr=0d0
 grad_elsr=0d0
!
 do i=1,npt ! charges
  __OUT(' Gridding charge #', i,' of ', npt);
  do img=-nimg,nimg ! go over charges in the original and in neighboring images
   xmin = x(i) + img*Lx - spt ; ! minimum location on grid for which potential nonzero
   xmax = xmin + spt + spt;
   imin = max (1, INT( (xmin-x0)*odx ) + 1); ! compute index limits on x-grid
   imax = min (nx,INT( (xmax-x0)*odx ) + 2);

! write(0,*) i, x(i), xmin, xmax, imin, imax
!stop
! proceed if charge supported
   if (imin<=imax) then
    xmin=xmin+spt;
    do jmg=-nimg,nimg
     ymin = y(i) + jmg*Ly - spt ; ! minimum location on grid for which potential nonzero
     ymax = ymin + spt + spt;
     jmin = max (1, INT( (ymin-y0)*ody ) + 1);
     jmax = min (ny,INT( (ymax-y0)*ody ) + 2);
! proceed if charge supported
     if (jmin<=jmax) then
      ymin=ymin+spt; ! center pt
      do kmg=-nimg,nimg
       zmin = z(i) + kmg*Lz - spt ; ! minimum location on grid for which potential nonzero
       zmax = zmin + spt + spt;
       kmin = max (1, INT( (zmin-z0)*odz ) + 1);
       kmax = min (nz,INT( (zmax-z0)*odz ) + 2);
! proceed if charge supported
       if (kmin<=kmax) then
        zmin=zmin+spt;
        do ii=imin,imax
         dx2=(xx(ii)-xmin)**2;
         do jj=jmin,jmax 
          dy2=dx2+(yy(jj)-ymin)**2;
          do kk=kmin,kmax
           r2=dy2+(zz(kk)-zmin)**2;
           if (r2>s2) cycle
           r=sqrt(r2)*oos; ! normalize by support
           __INCR(rho(ii,jj,kk), q(i) * fpoly(r))
          enddo !kk
         enddo !jj
        enddo !ii
       endif ! kmin
      enddo ! kmg
     endif ! jmin
    enddo ! jmg
   endif ! imin
  enddo !img
!
! also need to compute short-range potential energy at grid points
! do this in a primitive O(n^2) loop
  do j=i+1,npt
   do img=-nimg,nimg
    dx = (x(i)-x(j)) + Lx*img
    if (abs(dx)>spt) cycle
    do jmg=-nimg,nimg
     dy = (y(i)-y(j)) + Ly*jmg
     if (abs(dy)>spt) cycle
     dy2=dx**2 + dy**2
     if (dy2>s2) cycle
     do kmg=-nimg,nimg
      dz = (z(i)-z(j)) + Lz*kmg
      if (abs(dz)>spt) cycle
      dz2=dy2 + dz**2
      if (dz2 > s2 .or. dz2 < kzero) cycle
      r=sqrt(dz2)*oos
! short-range potential
      __INCR( elsr , q(i) * q(j) * fshort(r) ) ; ! might have contributions from multiple images

! short-range potential gradienst
      pre = q(i) * q(j) * fshortp(r) ;
!
      __INCR(elsr_dx(i) , pre*dx)
      __INCR(elsr_dy(i) , pre*dy)
      __INCR(elsr_dz(i) , pre*dz)
! reflect for other atom in pair
      __INCR(elsr_dx(j) , -pre*dx)
      __INCR(elsr_dy(j) , -pre*dy)
      __INCR(elsr_dz(j) , -pre*dz)
     enddo !kmg
    enddo !jmg
   enddo !img
  enddo ! j/npt
!
 enddo ! i/npt

! normalization
rho=rho/spt**3 ; ! filter support normalization
!
__SCALE(elsr, 1./(fourpi*eps)) ! no division by 2 because loop does not double count
__SCALE(elsr_dx, 1./(fourpi*eps))
__SCALE(elsr_dy, 1./(fourpi*eps))
__SCALE(elsr_dz, 1./(fourpi*eps))

end subroutine filt3p
