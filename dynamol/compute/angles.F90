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
      subroutine compute_angles3(e,alist,apar,r,fr,deriv) ! 3D
       use tlist
       use anglpar
       use constants
       implicit none
!
       type (toplist) :: alist
       type (angles) :: apar
       real*8 :: e
       real*8 :: r(:,:), fr(:,:)
       logical :: deriv
!
       real*8, parameter :: tol=1.0d-6
!
       real*8 :: theta, costh, sinth,&
& x1, x2, x3, y1, y2, y3,z1, z2, z3,&
& dx12, dy12, dz12, dx32, dy32, dz32,&
& vx, vy, vz, vn, ux, uy, uz, un,&
& fx12, fy12, fz12, fx32, fy32, fz32,&
& f1x, f1y, f1z, f2x, f2y, f2z, f3x, f3y, f3z, f,&
& fth, fub
!
       integer :: ind(4), i1, i2, i3, i4, i
       real*8 :: ktheta, theta0, kub, s0, dtheta
!do work:
       e=0d0
!
       do i=1, alist%last ! over all angles
!
        ind=alist%ind(:,i); i1=ind(1); i2=ind(2); i3=ind(3); i4=ind(4)
!
        x1=r(1,i1); x2=r(1,i2); x3=r(1,i3)
        y1=r(2,i1); y2=r(2,i2); y3=r(2,i3)
        z1=r(3,i1); z2=r(3,i2); z3=r(3,i3)
!
! deal with angle term first:
        dx12=x2-x1; dy12=y2-y1; dz12=z2-z1;
        dx32=x2-x3; dy32=y2-y3; dz32=z2-z3;
!
        ux=dx12; uy=dy12; uz=dz12;
        un=sqrt(ux*ux+uy*uy+uz*uz);
!
        vx=dx32; vy=dy32; vz=dz32;
        vn=sqrt(vx*vx+vy*vy+vz*vz);
!
        if (un.eq.0d0) un=1d0
        if (vn.eq.0d0) vn=1d0
!
        costh=(ux*vx+uy*vy+uz*vz)/(un*vn)
        sinth=sqrt(1d0-costh*costh)
        theta=atan2(sinth, costh) ! although acos would do just fine since we can't tell between theta & -theta
! add angle energy term
        ktheta=apar%ktheta(i4)
        theta0=apar%theta0(i4)
        kub=apar%kub(i4)
        s0=apar%s0(i4)
        dtheta=modulo(theta-theta0,twopi); if (dtheta.gt.pi) dtheta=dtheta-twopi; ! generic angle fix; wraps all angles
        fth=2.d0*ktheta*(dtheta) ! force prefactor
        e=e+ktheta*dtheta*dtheta
! add Urey-Bradley energy term
        if (kub.gt.0.) then
         dx12=x3-x1; dy12=y3-y1; dz12=z3-z1;
         theta=sqrt(dx12*dx12 + dy12*dy12 + dz12*dz12)
         dtheta=theta-s0
         fub=2.d0*kub*(dtheta) ! force prefactor
         e=e+kub*(dtheta)**2
        endif
!cccccccccccccccccc now compute derivative
        if (deriv) then
         un=1d0/un; ux=ux*un; uy=uy*un; uz=uz*un;
         vn=1d0/vn; vx=vx*vn; vy=vy*vn; vz=vz*vn;
!
         f=-fth/(max(sinth, tol)); ! avoid singularity at zero
!
         fx12=f*un*(vx-costh*ux);
         fy12=f*un*(vy-costh*uy);
         fz12=f*un*(vz-costh*uz);
!
         fx32=f*vn*(ux-costh*vx);
         fy32=f*vn*(uy-costh*vy);
         fz32=f*vn*(uz-costh*vz);
!
! compute derivatives w.r.t atom positions:
         f1x=-fx12; f1y=-fy12; f1z=-fz12;
         f2x= fx12+fx32; f2y= fy12+fy32; f2z= fz12+fz32;
         f3x=-fx32; f3y=-fy32; f3z=-fz32;
! add angle forces
         fr(1,i1)=fr(1,i1)-f1x; fr(2,i1)=fr(2,i1)-f1y; fr(3,i1)=fr(3,i1)-f1z
         fr(1,i2)=fr(1,i2)-f2x; fr(2,i2)=fr(2,i2)-f2y; fr(3,i2)=fr(3,i2)-f2z
         fr(1,i3)=fr(1,i3)-f3x; fr(2,i3)=fr(2,i3)-f3y; fr(3,i3)=fr(3,i3)-f3z
! add U/B forces
         if (kub.gt.0.) then
          if (theta.gt.tol) then
           f=fub/theta
          else
! nothing -- avoid singularity at near zero separation
          endif
          fx12=f*dx12; fy12=f*dy12; fz12=f*dz12
! compute derivatives w.r.t components
          f1x=-fx12; f1y=-fy12; f1z=-fz12;
          f2x= fx12; f2y= fy12; f2z= fz12;
! add angle forces
          fr(1,i1)=fr(1,i1)-f1x; fr(2,i1)=fr(2,i1)-f1y; fr(3,i1)=fr(3,i1)-f1z
          fr(1,i3)=fr(1,i3)-f2x; fr(2,i3)=fr(2,i3)-f2y; fr(3,i3)=fr(3,i3)-f2z
         endif ! U/B forces
        endif ! deriv
       enddo ! all angles
      end subroutine compute_angles3
