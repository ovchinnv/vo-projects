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
module sysmanip
 use output
 use system
 logical :: loud=.true.
 contains
  subroutine sysmanip_translate(dr,ind)
  implicit none
!
  character(len=19), parameter :: whoami='SYSMANIP_TRANSLATE'
!
  integer :: ind(:) ! indices of atoms (selection) whose properties to calculate; if first entry < 1, apply to all
  real*8 :: dr(ndim)
  integer :: i, n
  if (.not. system_coordinates_initialized) then
   call error(whoami, ' SYSTEM COORDINATES NOT INITIALIZED. ABORT.',0)
   return
  endif
!
  n=size(ind)
  if (n.eq.0) then ! empty selection
   call warning(whoami, ' EMPTY SELECTION.',0)
  else
   if ( ind(1).le.0 ) then ! use all indices
    do i=1,natom ; r(:,i)=r(:,i)+dr(:); enddo
   else ! a subset of indices
! are all indices valid?
    if (any(ind.lt.1).or.any(ind.gt.natom)) then
     call error(whoami, 'SOME INDICES IN SELECTION ARE INVALID. ABORT',0)
     return
    endif
    do i=1,n ; r(:,ind(i))=r(:,ind(i))+dr(:); enddo
   endif
  endif
!
  end subroutine sysmanip_translate
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine sysmanip_align_pc(qmass,ialign,imove)! align the principal components of the molecule with the Cartesian vectors
  use bestfit
  use constants
  use parser, only: ftoa, itoa
  implicit none
!
  logical :: qmass
  integer, dimension(:) :: ialign, imove
!
  character(len=17), parameter :: whoami='SYSMANIP_ALIGN_PC'
!
  integer :: i, j, n1, n2
  real*8 :: corr(3,3), rcom(3), evec(3,3), eval(3)
  real*8 :: mass(natom), rnew(3,natom)
  real*8 :: a, b, c, d, x, y, z, w
!
  if (ndim.lt.3) then
   call warning(whoami, ' ORIENTATION REQUIRES A 3D OBJECT. ABORT.',-1)
   return
  endif
!
  if (.not. system_coordinates_initialized) then
   call warning(whoami, ' SYSTEM COORDINATES NOT INITIALIZED. ABORT.',-1)
   return
  endif
!
  n1=size(ialign)
  if (n1.eq.0) then ! empty selection
   call warning(whoami, ' NO ATOMS SELECTED FOR ALIGNMENT.',-1)
   return
  endif
!
  n2=size(imove)
  if (n2.eq.0) then ! empty selection
   call warning(whoami, ' NO ATOMS SELECTED FOR ROTATION.',-1)
   return
  endif
!
  if (qmass) then ; mass=m ; else ; mass=1d0 ; endif ; d=sum(mass); if (d.gt.ERRTOL()) then ; d=1d0/d; else ; d=1d0 ; endif
  corr=0d0; rcom=0d0
  if (ialign(1).lt.0) then ! take all atom indices
   do i=1,natom
    x=r(1,i); y=r(2,i); z=r(3,i); w=mass(i)*d
!
    a=x*w; b=y*w; c=z*w;
    rcom(1)=rcom(1)+a; rcom(2)=rcom(2)+b; rcom(3)=rcom(3)+c;
!
    corr(1,1)=corr(1,1)+x*a; corr(1,2)=corr(1,2)+x*b; corr(1,3)=corr(1,3)+x*c;
                               corr(2,2)=corr(2,2)+y*b; corr(2,3)=corr(2,3)+y*c;
                                                           corr(3,3)=corr(3,3)+z*c;
   enddo
!
  else ! ialign
   do j=1,n1
    i=ialign(j)
    x=r(1,i); y=r(2,i); z=r(3,i); w=mass(i)*d
!
    a=x*w; b=y*w; c=z*w;
    rcom(1)=rcom(1)+a; rcom(2)=rcom(2)+b; rcom(3)=rcom(3)+c;
!
    corr(1,1)=corr(1,1)+x*a; corr(1,2)=corr(1,2)+x*b; corr(1,3)=corr(1,3)+x*c;
                               corr(2,2)=corr(2,2)+y*b; corr(2,3)=corr(2,3)+y*c;
                                                           corr(3,3)=corr(3,3)+z*c;
   enddo
!
  endif ! ialign
!
  do i=1,3; do j=i,3; corr(i,j)=corr(i,j)-rcom(i)*rcom(j) ; enddo; enddo;
  corr(2,1)=corr(1,2);
  corr(3,1)=corr(1,3);
  corr(3,2)=corr(2,3);
!
  call eig3s(corr,eval,evec)
!
  if (loud) then
   call message(whoami,'\t\t================== CENTER OF MASS ==================')
   call message(whoami,'   \t'//ftoa(rcom(1))//' '//ftoa(rcom(2))//' '//ftoa(rcom(3)))
   call message(whoami,'\t\t================= PRINCIPAL VECTORS ================')
   do i=1,3
    call message(whoami,'\t'//itoa(i)//': '//ftoa(evec(1,i))//' '//ftoa(evec(2,i))//' '//ftoa(evec(3,i)))
   enddo
   call message(whoami,'\t\t====================================================')
  endif
!
  if (imove(1).lt.0) then ! take all atom indices
   do i=1,3 ; r(i,:)=r(i,:)-rcom(i); enddo ! move to zero COM
   rnew=matmul(transpose(evec),r);
   do i=1,3 ; r(i,:)=rnew(i,:)+rcom(i); enddo ! move to original COM
  else
! use subindices: we are assuming that there are no repetitions
   do i=1,3 ; r(i,imove)=r(i,imove)-rcom(i); enddo ! move to zero COM
   r(:,imove)=matmul(transpose(evec),r(:,imove));
   do i=1,3 ; r(i,imove)=r(i,imove)+rcom(i); enddo ! move to zero COM
  endif
!
  end subroutine sysmanip_align_pc
!
end module sysmanip
