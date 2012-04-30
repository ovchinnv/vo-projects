     interface
      subroutine compute_bonds3(e,blist,bpar,r,fr,deriv) ! 3D
       use tlist
       use bondpar
       implicit none
!
       type (toplist) :: blist
       type (bonds) :: bpar
       float :: e
       float :: r(:,:), fr(:,:)
       bool :: deriv
      end subroutine compute_bonds3
!
      subroutine compute_angles3(e,alist,apar,r,fr,deriv) ! 3D
       use tlist
       use anglpar
       implicit none
!
       type (toplist) :: alist
       type (angles) :: apar
       float :: e
       float :: r(:,:), fr(:,:)
       bool :: deriv
      end subroutine compute_angles3
!
      subroutine compute_dihes3(e,dlist,dpar,r,fr,deriv) ! 3D
       use tlist
       use dihepar
       use constants, only: pi, twopi
       implicit none
!
       type (toplist) :: dlist
       type (dihes) :: dpar
       float :: e
       float :: r(:,:), fr(:,:)
       bool :: deriv
      end subroutine compute_dihes3
     end interface
