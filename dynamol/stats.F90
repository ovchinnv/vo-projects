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
! subroutines to compute various statistical quantities
module stats
 use constants
 contains
 function calcKE(v, m)
 implicit none
 integer :: natom, ndim
 real*8 :: v(:,:), m(:)
 real*8, dimension (2) :: calcKE
! assume that the dimensions are correct for multiplication to make sense
 natom=size(v,2)
 ndim=size(v,1)
!
 calcKE(1)=0.5d0*sum(matmul(v**2,m))
 calcKE(2)=calcKE(1)*2/ndim/kboltzmann/natom
!
 end function calcKE
end module stats
