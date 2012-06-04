/*#define __WRN(__WHO,__MSG) write(0,*) 'WARNING FROM: ',__WHO,': ',__MSG*/
/*#define __PRINT(__MSG) write(0,'(A)') __MSG*/
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
      subroutine compute_dihes3(e,dlist,dpar,r,fr,deriv) ! 3D
       use tlist
       use dihepar
       use constants
       implicit none
!
       type (toplist) :: dlist
       type (dihes) :: dpar
       real*8 :: e
       real*8 :: r(:,:), fr(:,:)
       logical :: deriv
!
       real*8, parameter :: tol=1.0d-6
!
       real*8 :: theta, costh, sinth,&
& x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4,&
& dx12, dy12, dz12, dx23, dy23, dz23, dx34, dy34, dz34,&
& vx, vy, vz, vn, ux, uy, uz, un, wx, wy, wz, wn,&
& dcosdux, dcosduy, dcosduz, dcosdvx, dcosdvy, dcosdvz,&
& dsindux, dsinduy, dsinduz, dsindwx, dsindwy, dsindwz,&
& fx12, fy12, fz12, fx23, fy23, fz23, fx34, fy34, fz34,&
& f1x, f1y, f1z, f2x, f2y, f2z, f3x, f3y, f3z,&
& f4x, f4y, f4z, f
!
       integer :: ind(5), i1, i2, i3, i4, i5, i, mult
       real*8 :: kchi, dtheta, delta
!do work:
       e=0d0
!
       do i=1, dlist%last ! over all angles
!
        ind=dlist%ind(:,i); i1=ind(1); i2=ind(2); i3=ind(3); i4=ind(4); i5=ind(5)
!
        x1=r(1,i1); x2=r(1,i2); x3=r(1,i3); x4=r(1,i4)
        y1=r(2,i1); y2=r(2,i2); y3=r(2,i3); y4=r(2,i4)
        z1=r(3,i1); z2=r(3,i2); z3=r(3,i3); z4=r(3,i4)
!
        dx12=x2-x1; dy12=y2-y1; dz12=z2-z1;
        dx23=x3-x2; dy23=y3-y2; dz23=z3-z2;
        dx34=x4-x3; dy34=y4-y3; dz34=z4-z3;
! note:
! costh = [ (d34 x d32) . (d23 x d21) / |(d34 x d23)| |(d23 x d12)| ] =
! = acos [ (u . v) / |u| |v| ] (definition of u and v)
        ux=dy23*dz34-dy34*dz23;
        uy=dz23*dx34-dz34*dx23;
        uz=dx23*dy34-dx34*dy23;
        un=sqrt(ux*ux+uy*uy+uz*uz);
!
        vx=dy12*dz23-dy23*dz12;
        vy=dz12*dx23-dz23*dx12;
        vz=dx12*dy23-dx23*dy12;
        vn=sqrt(vx*vx+vy*vy+vz*vz);
!
        wx=dy23*vz-vy*dz23;
        wy=dz23*vx-vz*dx23;
        wz=dx23*vy-vx*dy23;
        wn=sqrt(wx*wx+wy*wy+wz*wz);
!
        if (un.eq.0d0) un=1d0
        if (vn.eq.0d0) vn=1d0
        if (wn.eq.0d0) wn=1d0
        costh=(ux*vx+uy*vy+uz*vz)/(un*vn)
        sinth=(wx*ux+wy*uy+wz*uz)/(wn*un)
        theta=atan2(sinth, costh)
!
! add dihedral energy term
!V(dihedral) = Kchi(1 + cos(n(chi) - delta))
!
!Kchi: kcal/mole
!n: multiplicity
!delta: degrees
        kchi=dpar%kchi(i5)
        delta=dpar%delta(i5)
        mult=dpar%mult(i5)
        if (mult.gt.0) then
         theta=mult*theta
         dtheta=modulo(theta-delta,twopi); if (dtheta.gt.pi) dtheta=dtheta-twopi; ! generic angle fix; wraps all angles
         f=-mult*kchi*sin(dtheta) ! force prefactor
         e=e+kchi*(1+cos(dtheta))
        else ! 0 multiplicity: assume this is an improper dihedral
         dtheta=modulo(theta-delta,twopi); if (dtheta.gt.pi) dtheta=dtheta-twopi; ! generic angle fix; wraps all angles
         f=2.0*kchi*dtheta ! force prefactor
         e=e+kchi*dtheta*dtheta
        endif
!cccccccccccccccccc now compute derivative
        if (deriv) then
         un=1d0/un
         ux=ux*un; uy=uy*un; uz=uz*un;
! decide which form to use: 1/sin(theta) is singular for theta=n*pi
         if (abs(sinth).gt.0.1) then ! use simpler version (1/sin defined)
          vn=1d0/vn
          vx=vx*vn; vy=vy*vn; vz=vz*vn;
          dcosdux=un*(vx-costh*ux);
          dcosduy=un*(vy-costh*uy);
          dcosduz=un*(vz-costh*uz);
!
          dcosdvx=vn*(ux-costh*vx);
          dcosdvy=vn*(uy-costh*vy);
          dcosdvz=vn*(uz-costh*vz);
! calculate derivative w.r.t. `bonds' 12 23 34;
          f=-f/sinth;
          fx12 = f*(dy23*dcosdvz-dz23*dcosdvy);
          fy12 = f*(dz23*dcosdvx-dx23*dcosdvz);
          fz12 = f*(dx23*dcosdvy-dy23*dcosdvx);
!
          fx23 = f*(dy34*dcosduz-dz34*dcosduy- (dy12*dcosdvz-dz12*dcosdvy) )
          fy23 = f*(dz34*dcosdux-dx34*dcosduz- (dz12*dcosdvx-dx12*dcosdvz) )
          fz23 = f*(dx34*dcosduy-dy34*dcosdux- (dx12*dcosdvy-dy12*dcosdvx) )
!
          fx34 = f*(-dy23*dcosduz+dz23*dcosduy);
          fy34 = f*(-dz23*dcosdux+dx23*dcosduz);
          fz34 = f*(-dx23*dcosduy+dy23*dcosdux);
!
!
         else ! use less efficient form (1/cos defined)
          wn=1d0/wn
          wx=wx*wn; wy=wy*wn; wz=wz*wn;
!
          dsindwx=wn*(ux-sinth*wx);
          dsindwy=wn*(uy-sinth*wy);
          dsindwz=wn*(uz-sinth*wz);
!
          dsindux=un*(wx-sinth*ux);
          dsinduy=un*(wy-sinth*uy);
          dsinduz=un*(wz-sinth*uz);
!
          f=f/costh;
          fx12=f*(dsindwx*(dy23*dy23+dz23*dz23)&
& -dsindwy*(dx23*dy23)&
& -dsindwz*(dx23*dz23))
          fy12=f*(-dsindwx*(dx23*dy23)&
& + dsindwy*(dx23*dx23+dz23*dz23)&
& - dsindwz*(dy23*dz23))
          fz12=f*(-dsindwx*(dx23*dz23)&
& - dsindwy*(dy23*dz23)&
& + dsindwz*(dy23*dy23+dx23*dx23))
!
          fx23=f*(dy34*dsinduz-dz34*dsinduy &
& -dsindwx*(dz23*dz12+dy23*dy12) &
& +dsindwy*(2d0*dy12*dx23-dy23*dx12)&
& +dsindwz*(2d0*dz12*dx23-dz23*dx12))
          fy23=f*(dz34*dsindux-dx34*dsinduz &
& +dsindwx*(2d0*dx12*dy23-dx23*dy12) &
& -dsindwy*(dx12*dx23+dz12*dz23) &
& +dsindwz*(2d0*dz12*dy23-dz23*dy12))
          fz23=f*(dx34*dsinduy-dy34*dsindux &
& +dsindwx*(2d0*dx12*dz23-dx23*dz12)&
& +dsindwy*(2d0*dy12*dz23-dz12*dy23)&
& -dsindwz*(dx12*dx23+dy12*dy23))
!
          fx34=f*(dsinduy*dz23-dsinduz*dy23)
          fy34=f*(dsinduz*dx23-dsindux*dz23)
          fz34=f*(dsindux*dy23-dsinduy*dx23)
         endif ! sinth
! compute derivatives w.r.t positions:
         f1x=-fx12; f1y=-fy12; f1z=-fz12;
         f2x= fx12-fx23; f2y= fy12-fy23; f2z= fz12-fz23;
         f3x= fx23-fx34; f3y= fy23-fy34; f3z= fz23-fz34;
         f4x= fx34; f4y= fy34; f4z= fz34;
! add dihedral forces
         fr(1,i1)=fr(1,i1)-f1x; fr(2,i1)=fr(2,i1)-f1y; fr(3,i1)=fr(3,i1)-f1z
         fr(1,i2)=fr(1,i2)-f2x; fr(2,i2)=fr(2,i2)-f2y; fr(3,i2)=fr(3,i2)-f2z
         fr(1,i3)=fr(1,i3)-f3x; fr(2,i3)=fr(2,i3)-f3y; fr(3,i3)=fr(3,i3)-f3z
         fr(1,i4)=fr(1,i4)-f4x; fr(2,i4)=fr(2,i4)-f4y; fr(3,i4)=fr(3,i4)-f4z
!
        endif ! deriv
       enddo ! all angles
      end subroutine compute_dihes3
