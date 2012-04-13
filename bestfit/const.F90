! This source file was was generated automatically from a master source
! code tree, which may or may not be distributed with this code.
! If you edit this file and (rather than the master source file)
! your changes will be lost after the next pull from the master tree.
! This approach makes it possible to have the same master source code
! interfaced with several applications at the source level
!
module const
  private
  !******************** Victor Ovchinnikov 1/2012: rough estimate of machine precision
  public ERRTOL
  real*8, save, private :: ERRTOL_ = -1 ! error tolerance
 contains
  function ERRTOL()
   implicit none
   real*8 :: a, b, ERRTOL
   if (ERRTOL_.lt.0) then
    a=1.0; b=1.0;
    do
     if ( a - a / b .ne. a ) then ; b = b * 10 ; else ; ERRTOL_ = 50 / b ; exit ; endif
    enddo
   endif
   ERRTOL=ERRTOL_
  end function ERRTOL
  !
  !******************** End Victor Ovchinnikov 1/2012: rough estimate of machine precision
  !
end module const
