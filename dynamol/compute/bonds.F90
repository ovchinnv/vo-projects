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
      subroutine compute_bonds3(e,blist,bpar,r,fr,deriv) ! 3D
       use tlist
       use bondpar
       implicit none
!
       type (toplist) :: blist
       type (bonds) :: bpar
       real*8 :: e
       real*8 :: r(:,:), fr(:,:)
       logical :: deriv
!
       real*8, parameter :: tol=1.0d-6
!
       real*8 :: theta,&
& x1, x2, y1, y2, z1, z2,&
& dx12, dy12, dz12,&
& fx12, fy12, fz12,&
& f1x, f1y, f1z, f2x, f2y, f2z, f
!
       integer :: ind(3), i1, i2, i3, i
       real*8 :: kb, b0
!do work:
       e=0d0
!
       do i=1, blist%last ! over all angles
!
        ind=blist%ind(:,i); i1=ind(1); i2=ind(2); i3=ind(3);
!
        x1=r(1,i1); x2=r(1,i2);
        y1=r(2,i1); y2=r(2,i2);
        z1=r(3,i1); z2=r(3,i2);
!
        dx12=x2-x1; dy12=y2-y1; dz12=z2-z1;
        theta=sqrt(dx12*dx12 + dy12*dy12 + dz12*dz12)
! add energy term
        kb=bpar%kb(i3)
        b0=bpar%b0(i3)
        f=2.d0*kb*(theta-b0) ! force prefactor
        e=e+kb*(theta-b0)**2
!cccccccccccccccccc now compute derivative
        if (deriv) then
         if (theta.gt.tol) then
          f=f/theta
         else
! nothing : avoid singularity at near zero separation
         endif
         fx12=f*dx12; fy12=f*dy12; fz12=f*dz12
! compute derivatives w.r.t components
         f1x=-fx12; f1y=-fy12; f1z=-fz12;
         f2x= fx12; f2y= fy12; f2z= fz12;
! add angle forces
         fr(1,i1)=fr(1,i1)-f1x; fr(2,i1)=fr(2,i1)-f1y; fr(3,i1)=fr(3,i1)-f1z
         fr(1,i2)=fr(1,i2)-f2x; fr(2,i2)=fr(2,i2)-f2y; fr(3,i2)=fr(3,i2)-f2z
!
        endif ! deriv
       enddo ! all angles
      end subroutine compute_bonds3
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
