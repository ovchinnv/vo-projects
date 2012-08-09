/*#define __WRN(__WHO,__MSG) write(0,*) 'WARNING FROM: ',__WHO,': ',__MSG*/
/*#define __PRINT(__MSG) write(0,'(A)') __MSG*/
/*#define __PRINT(__MSG) call plainmessage(__MSG)*/
/*#define __PRINTL(__MSG,__LEVEL) call plainmessage(__MSG,__LEVEL)*/
/*COORDINATES AND MASSES:*/
/*#define __INDX(__STR, __STRLEN, __TEST, __TESTLEN)  index(__STR(1:min(__STRLEN,len(__STR))),__TEST(1:min(__TESTLEN,len(__TEST))))*/
/*
#ifdef __IMPNONE
#undef __IMPNONE
#endif
#define __IMPNONE
*/
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
!
! CV_DIST_COM.MOD
!
! ROUTINES FOR COLLECTIVE VARIABLE `DISTANCE' BETWEEN CENTERS OF MASS OF ATOMS
!
!**CHARMM_ONLY**!##IF STRINGM
!
      module cv_dist_com
      use cv_common ! `base' cv class
      use cv_types
      use ivector ! vector class for storing atom lists
!
      ! subroutines
      public cv_dist_com_add
      public cv_dist_com_calc
      public cv_dist_com_list
!
      contains
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       function cv_dist_com_add(atom_list,k,gamma,weight) ! note: i is the atom index in the PSF
       use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
       implicit none
       real*8 :: k, gamma, weight
       type (int_vector), dimension(2) :: atom_list
! locals
       integer :: i, j, l, m, ind, num_int, ncom(2)
       logical :: found, cv_dist_com_add
       character(len=17) :: whoami
       type (int_vector) :: unique_amap_ptr ! contains a unique list of atom map pointers (to speed up gradient computation, etc)
!c
       data whoami/' CV_DIST_COM_ADD>'/
!
! check for duplicate CV (exact identical entry only)
       found=.false.
       do l=1, cv%num_cv
        if (cv%type(l).eq.dist_com) then
         found=.true.
         do i=1,2
          ncom(i)=atom_list(i)%last
          if (found) found=(ncom(i).eq.cv%priv(l)%p(i))
         enddo
         ind=3
         do i=1,2
          do j=1,ncom(i)
           if (found) found= &
     & (atom_list(i)%i(j).eq.cv%amap%i(cv%priv(l)%p(ind)))
           ind=ind+1
          enddo
         enddo
        endif
        if (found) exit
       enddo
!
       if (.not.found) then ! (if found -- do nothing)
        l=cv_common_add(k,gamma,weight,dist_com) ! get a new cv index
        if (l.gt.0) then
! allocate private data
! space needed:
         do i=1,2 ; ncom(i)=atom_list(i)%last ; enddo
         num_int = 2 + sum(ncom) ! number of ints needded for storage
!
         allocate(cv%priv(l)%p(num_int));
         cv%priv(l)%p(1:2)=ncom(1:2)
! now add atom indices
         ind=3
         do i=1,2
          do j=1,ncom(i)
           m=atom_list(i)%i(j)
           if (m.le.0) call warning(whoami, ' INVALID ATOM INDEX SPECIFIED.', 0)
           cv%priv(l)%p(ind)=int_vlist_uaddu(cv%amap,m,l) ! add indices into unique map
           m=int_vector_uadd(unique_amap_ptr,cv%priv(l)%p(ind))
           ind=ind+1
          enddo
         enddo
!
         m=unique_amap_ptr%last ! number of unique atoms this cv depends on
         allocate(cv%priv(l)%amap_ptr(m+1)) ! add one to include length of list
         cv%priv(l)%amap_ptr(1)=m;
         cv%priv(l)%amap_ptr(2:m+1)=unique_amap_ptr%i(1:m)
         call int_vector_done(unique_amap_ptr)
!
         cv_dist_com_add=.true.
        else ! out of bounds
         call warning(whoami, ' ERROR ADDING DIST_COM CV. NOTHING DONE.', 0)
         cv_dist_com_add=.false.
        endif
       else ! found
         call warning(whoami, ' DIST_COM CV ALREADY PRESENT. NOTHING DONE.', 0)
         cv_dist_com_add=.false.
       endif
       end function cv_dist_com_add
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine cv_dist_com_calc(i,x,y,z,mass,fx,fy,fz, &
     & calctheta,deriv,addforce,fext)
       implicit none
!
       real*8 :: x(:), y(:), z(:), &
     & fx(:), fy(:), fz(:), mass(:)
       real*8, optional :: fext ! `external' force for planar dynamics (see smcv_addforce in smcv_master)
       integer :: i ! cv index
       logical :: calctheta ! whether or not to calculate theta(x); if not, use theta=cv%r(i,instant)
                            ! note that if calctheta=.false., we don't calculate the derivatives!
       logical :: deriv ! whether or not to calculate derivatives w.r.t. x
       logical :: addforce ! whether or not to add forces on simulation atoms
!
       real*8, parameter :: tol=1.0e-10
       real*8 :: dummy
!
       integer, allocatable, dimension(:) :: ind1, ind2
! & ,indpsf1, indpsf2
       real*8, allocatable, dimension(:) :: x1, x2
       real*8, allocatable, dimension(:) :: y1, y2
       real*8, allocatable, dimension(:) :: z1, z2
       real*8, allocatable, dimension(:) :: m1, m2
!
       integer :: ncom1, ncom2, ind ! , indpsf
       integer :: j, ii, jj
! variables for cv and derivative calculations
       real*8 :: theta, & ! the convention is that any CV value is called `theta'
     & xcom1, xcom2, ycom1, ycom2, zcom1, zcom2, &
     & dx12, dy12, dz12, fx12, fy12, fz12, &
     & f1x, f1y, f1z, f2x, f2y, f2z, f, w, totm1, totm2
       integer, pointer :: priv(:)
! do work:
! look up atom index in the private array of CV i; obtain the PSF index from amap
       ncom1=cv%priv(i)%p(1)
       ncom2=cv%priv(i)%p(2)
!
       allocate(ind1(ncom1),x1(ncom1),y1(ncom1),z1(ncom1),m1(ncom1))
       allocate(ind2(ncom2),x2(ncom2),y2(ncom2),z2(ncom2),m2(ncom2))
! allocate(indpsf1(ncom1),indpsf2(ncom2))
!
! extract indices into the atom map
       ii=3; jj=ii+ncom1-1; ind1=cv%priv(i)%p(ii:jj)
       ii=jj+1; jj=ii+ncom2-1; ind2=cv%priv(i)%p(ii:jj)
!
       do j=1, ncom1;
        ind=cv%amap%i(ind1(j)) ! actual psf index
! indpsf1(j)=ind ;
                      x1(j)=x(ind) ; y1(j)=y(ind) ; z1(j)=z(ind) ;
                      m1(j)=mass(ind)
       enddo
       do j=1, ncom2;
        ind=cv%amap%i(ind2(j))
! indpsf2(j)=ind;
                      x2(j)=x(ind); y2(j)=y(ind) ; z2(j)=z(ind)
                      m2(j)=mass(ind)
       enddo
! now we have all the coordinates needed to compute distance and its derivative :
       if (calctheta) then
!
        totm1=1d0/sum(m1); totm2=1d0/sum(m2);
! normalize masses
        m1=m1*totm1; m2=m2*totm2;
!
        xcom1=0.; xcom2=0.;
        ycom1=0.; ycom2=0.;
        zcom1=0.; zcom2=0.;
        do j=1, ncom1;
         xcom1=xcom1+x1(j)*m1(j);
         ycom1=ycom1+y1(j)*m1(j);
         zcom1=zcom1+z1(j)*m1(j);
        enddo
        do j=1, ncom2;
         xcom2=xcom2+x2(j)*m2(j);
         ycom2=ycom2+y2(j)*m2(j);
         zcom2=zcom2+z2(j)*m2(j);
        enddo
! compute centers of mass
! xcom1=xcom1*totm1; ycom1=ycom1*totm1; zcom1=zcom1*totm1;
! xcom2=xcom2*totm2; ycom2=ycom2*totm2; zcom2=zcom2*totm2;
!
        dx12=xcom2-xcom1; dy12=ycom2-ycom1; dz12=zcom2-zcom1;
!
        theta=sqrt(dx12*dx12 + dy12*dy12 + dz12*dz12)
        cv%r(i,previnst)=cv%r(i,instant)
        cv%r(i,instant)=theta
!cccccccccccccccccc now compute derivative
        if (deriv) then
         if (theta.gt.tol) then
          f=1d0/theta
         else
          f=1d0 ! avoid singularity at near zero separation
         endif
         fx12=f*dx12; fy12=f*dy12; fz12=f*dz12
! compute derivatives w.r.t COMs:
         f1x=-fx12; f1y=-fy12; f1z=-fz12;
         f2x= fx12; f2y= fy12; f2z= fz12;
! compute derivatives w.r.t. COM components
! initialize derivative arrays
! do j=1,cv%amap%last
         priv=>cv%priv(i)%amap_ptr ! this should copy the pointers statically; the pointers point to the _same_ data
         do jj=2,priv(1)+1 ! only a subset of indices needs to be considered
          j=priv(jj)
!
          do ii=1,2
          cv%gradx(i,j,ii)=0d0;cv%grady(i,j,ii)=0d0;cv%gradz(i,j,ii)=0d0
          enddo
         enddo
! loop over all groups
! COM1
! loop over all relevant indices in the map
         do j=1,ncom1
          ind=ind1(j)
          w=sqrt(totm1/m1(j)) ! for mass-weighting
!
          f=f1x*m1(j)
          cv%gradx(i,ind,1)=f ! derivative w.r.t x component of atom indpsf1(j)
          cv%gradx(i,ind,2)=f*w
!
          f=f1y*m1(j)
          cv%grady(i,ind,1)=f
          cv%grady(i,ind,2)=f*w
!
          f=f1z*m1(j)
          cv%gradz(i,ind,1)=f
          cv%gradz(i,ind,2)=f*w
         enddo
! COM2
         do j=1,ncom2
          ind=ind2(j)
          w=sqrt(totm2/m2(j))
!
          f=f2x*m2(j)
          cv%gradx(i,ind,1)=cv%gradx(i,ind,1)+f
          cv%gradx(i,ind,2)=cv%gradx(i,ind,2)+f*w
!
          f=f2y*m2(j)
          cv%grady(i,ind,1)=cv%grady(i,ind,1)+f
          cv%grady(i,ind,2)=cv%grady(i,ind,2)+f*w
!
          f=f2z*m2(j)
          cv%gradz(i,ind,1)=cv%gradz(i,ind,1)+f
          cv%gradz(i,ind,2)=cv%gradz(i,ind,2)+f*w
         enddo
        endif ! deriv
       else ! .not.calctheta
        theta=cv%r(i,instant)
       endif
!
! NOTE that the forces calculated here are NOT acting on the CV, as in
! the evolution subroutine, but on the simulation atoms
       dummy=cv%r(i,zcur)-theta ! zcur contains reference coords (combination of main+comp)
!
       f=cv%k(i)*dummy ! note that f is not divided by 1/sqrt(mass)
       cv%r(i,forces2)=f
       if (addforce) then
        if (present(fext)) f=fext ! override for planar sampling
! loop over all groups
! COM1
! loop over all relevant indices in the map
! note: fx are actually components of grad V, not -grad V!
! do j=1,ncom1
! ind=ind1(j)
! indpsf=indpsf1(j) ! index into psf
! fx(indpsf)=fx(indpsf)-f*cv%gradx(i,ind,1)
! fy(indpsf)=fy(indpsf)-f*cv%grady(i,ind,1)
! fz(indpsf)=fz(indpsf)-f*cv%gradz(i,ind,1)
! enddo
! COM2
! do j=1,ncom2
! ind=ind2(j)
! indpsf=indpsf2(j) ! index into psf
! fx(indpsf)=fx(indpsf)-f*cv%gradx(i,ind,1)
! fy(indpsf)=fy(indpsf)-f*cv%grady(i,ind,1)
! fz(indpsf)=fz(indpsf)-f*cv%gradz(i,ind,1)
! enddo
! more efficient and conceptually correct way:
        priv=>cv%priv(i)%amap_ptr
        do ii=2, priv(1)+1
         jj=priv(ii)
         j=cv%amap%i(jj) ! psf index
         fx(j)=fx(j)-f*cv%gradx(i,jj,1)
         fy(j)=fy(j)-f*cv%grady(i,jj,1)
         fz(j)=fz(j)-f*cv%gradz(i,jj,1)
        enddo
       endif ! addforce
! free memory
       deallocate(ind1, ind2, & ! indpsf1, indpsf2,
     & x1, x2, y1, y2, z1, z2,m1, m2)
       end subroutine cv_dist_com_calc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine cv_dist_com_list(i)
       use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
       use multicom_aux
       use mpi
       use psf, only : atoms
       implicit none
 character(len=200) :: msg___(10)=(/'','','','','','','','','',''/); integer :: i_
!
       integer :: i, j, type, ii, jj, iatom
       character(len=8) :: sid, rid, ren, ac
       character(len=18) :: whoami
!
       integer :: ncom1, ncom2
       integer, allocatable, dimension(:) :: ind1, ind2
       data whoami /' CV_DIST_COM_LIST>'/
!
!
       if (MPI_COMM_STRNG.eq.MPI_COMM_NULL) return ! only replica heads stay
!
! check type just in case
       type=cv%type(i)
       if (type.ne.dist_com) then
        call warning(whoami, ' WRONG CV TYPE RECEIVED.', 0)
       endif
!
       if (ME_STRNG.eq.0) then
        ncom1=cv%priv(i)%p(1)
        ncom2=cv%priv(i)%p(2)
!
        allocate(ind1(ncom1),ind2(ncom2))
!
! extract indices into the atom map
        ii=3; jj=ii+ncom1-1; ind1=cv%priv(i)%p(ii:jj)
        ii=jj+1; jj=ii+ncom2-1; ind2=cv%priv(i)%p(ii:jj)
!
        write(msg___,'(A)') '\t DISTANCE-COM, GROUP 1' ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
        do j=1, ncom1;
         iatom=cv%amap%i(ind1(j)) ! actual psf index
         sid=atoms%segid(iatom); rid=atoms%resid(iatom); ren=atoms%resname(iatom); ac=atoms%aname(iatom);
         write(msg___,667) '\t',j, iatom, sid, rid, ren, ac ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
        enddo
        write(msg___,'(A)') '\t DISTANCE-COM, GROUP 2' ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
        do j=1, ncom2;
         iatom=cv%amap%i(ind2(j))
         sid=atoms%segid(iatom); rid=atoms%resid(iatom); ren=atoms%resname(iatom); ac=atoms%aname(iatom);
         write(msg___,667) '\t',j, iatom, sid, rid, ren, ac ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
        enddo
        deallocate(ind1,ind2)
       endif
!
 667 format(A,2I8,' ',4A)
!
       end subroutine cv_dist_com_list
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end module cv_dist_com
!
!**CHARMM_ONLY**!##ENDIF
