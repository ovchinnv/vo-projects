# 1 "const.ftn"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "const.ftn"
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
