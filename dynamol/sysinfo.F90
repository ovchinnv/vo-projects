! **********************************************************************!
! This source file was was generated automatically from a master source !
! code tree, which may or may not be distributed with this code, !
! because it is up to the distributor, and not up to me. !
! If you edit this file (rather than the master source file) !
! your changes will be lost if another pull from the master tree occurs.!
! In case you are wondering why, this approach makes it possible for !
! me to have the same master source code interfaced with different !
! applications (some of which are written in a way that is quite far !
! from being object-oriented) at the source level. !
! **********************************************************************!
module sysinfo
 use output
 use system
! routines to obtain various properties of the simulation system belong here
 public sysinfo_dimens ! return molecule dimensions: coordinates of center, minimum and maximum coordinate values
!
 contains
!***************************************************************************************************
  function sysinfo_dimens(qmass,ind)
  implicit none
!
  character(len=14), parameter :: whoami='SYSINFO_DIMENS'
!
  integer :: ind(:) ! indices of atoms (selection) whose properties to calculate; if missing do all
  logical :: qmass ! whether to use mass weighting
  real*8 :: dimens(ndim,3)=0, sysinfo_dimens(ndim,3)
  real*8 :: oomtot
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
!
    if (qmass) then
     oomtot=1d0/sum(m(1:natom));
     do i=1,ndim ; dimens(i,1)=dot_product( r(i,1:natom) , m(1:natom) )*oomtot ; enddo; ! COM
    else
     do i=1,ndim ; dimens(i,1)=sum( r(i,1:natom) )/natom ; enddo ! COG
    endif
!
    do i=1,ndim ;
     dimens(i,2)=minval(r(i,1:natom));
     dimens(i,3)=maxval(r(i,1:natom));
    enddo
!
   else ! a subset of indices
! are all indices valid?
    if (any(ind.le.0).or.any(ind.gt.natom)) then
     call error(whoami, 'SOME INDICES IN SELECTION ARE INVALID. ABORT',0)
     return
    endif
!
    if (qmass) then
     oomtot=1d0/sum(m(ind(:)));
     do i=1,ndim ; dimens(i,1)=dot_product( r(i,ind) , m(ind) )*oomtot ; enddo; ! COM
    else
     do i=1,ndim ; dimens(i,1)=sum( r(i,ind) )/natom ; enddo ! COG
    endif
!
    do i=1,ndim ;
     dimens(i,2)=minval(r(i,ind));
     dimens(i,3)=maxval(r(i,ind));
    enddo
   endif
!
  endif
!
  sysinfo_dimens=dimens
!
  end function sysinfo_dimens
!
end module sysinfo
