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
      module bondpar
! define a derived type to store bond parameters
      type bonds
       character(len=8), dimension(:), pointer :: a1, a2
       real*8, dimension(:), pointer :: kb, b0
! V(bond) = Kb(b - b0)**2
!Kb: kcal/mole/A**2
!b0: A
       integer :: length ! length of the vector
       integer :: last ! index of last element
       logical :: initialized=.false. ! has the vector been initialized
      end type bonds
!
      integer, parameter, private :: expand_incr=200
!
      contains
!******************************************** bond data routines*********************
       subroutine bonds_init( v )
       implicit none
       type (bonds) :: v
! if (associated(v%...)) deallocate(v%...) ! testing unassigned pointer is an error!
       allocate(v%a1(expand_incr),v%a2(expand_incr),v%kb(expand_incr),v%b0(expand_incr))
!
       v%a1=''; v%a2=''; v%kb=0; v%b0=0
       v%length=expand_incr
       v%last=0
       v%initialized=.true.
       end subroutine bonds_init
!ccccc
       subroutine bonds_done( v )
       implicit none
       type (bonds) :: v
       if (associated(v%a1)) deallocate(v%a1)
       if (associated(v%a2)) deallocate(v%a2)
       if (associated(v%kb)) deallocate(v%kb)
       if (associated(v%b0)) deallocate(v%b0)
       v%length=0
       v%last=0
       v%initialized=.false.
       end subroutine bonds_done
!
       subroutine bonds_expand( v )
       implicit none
       type (bonds) :: v
       integer :: newlength
       character(len=8), dimension(:), pointer :: a1, a2
       real*8, dimension(:), pointer :: kb, b0
!
       if (.not.v%initialized) then
        call bonds_init(v)
       else
! assume length is valid
        newlength=v%length+expand_incr ! temporary storage space
        allocate(a1(newlength),a2(newlength),kb(newlength),b0(newlength)) ! copy old data
        a1(1:v%length)=v%a1; a2(1:v%length)=v%a2; kb(1:v%length)=v%kb; b0(1:v%length)=v%b0
        deallocate(v%a1,v%a2,v%kb,v%b0)
        v%a1=>a1; v%a2=>a2; v%kb=>kb; v%b0=>b0
        v%length=newlength
       endif
       end subroutine bonds_expand
!ccccc
       function bonds_add(v,a1,a2,kb,b0)
       implicit none
       type (bonds) :: v
       integer :: bonds_add
       real*8 :: kb, b0
       character(len=*) :: a1, a2
       integer :: j
!
       if (.not.v%initialized) call bonds_init(v)
! add element to the list
       if (v%last.eq.v%length) call bonds_expand(v)
       j=v%last+1
       v%a1(j)=a1; v%kb(j)=kb
       v%a2(j)=a2; v%b0(j)=b0
       v%last=j
       bonds_add=j
       end function bonds_add
!ccccc
       function bonds_uadd(v,a1,a2,kb,b0)
       use output, only: warning
       implicit none
       type (bonds) :: v
       integer :: bonds_uadd
       real*8 :: kb, b0
       character(len=*) :: a1, a2
       integer :: j
!
       if (.not.v%initialized) call bonds_init(v)
       do j=1,v%last
        if ( (v%a1(j).eq.a1.and.v%a2(j).eq.a2).or.(v%a1(j).eq.a2.and.v%a2(j).eq.a1) ) then
! found element
         bonds_uadd=j
         call warning('BONDS_UADD','Parameters for bond "'//v%a1(j)//'--'//v%a2(j)//'" already present. Will overwrite.',0)
         v%kb(j)=kb; v%b0(j)=b0
         return
        endif
       enddo
! add element to the list: use regular routine
       bonds_uadd=bonds_add(v,a1,a2,kb,b0)
!
       end function bonds_uadd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       function bonds_getind( v,a1,a2 )
       implicit none
       type (bonds) :: v
       integer :: bonds_getind
       character(len=*) :: a1, a2
       integer :: j
!
       bonds_getind=-1
       if (.not.v%initialized) then
!
       else
        do j=1,v%last
         if ( (v%a1(j).eq.a1.and.v%a2(j).eq.a2).or.(v%a1(j).eq.a2.and.v%a2(j).eq.a1) ) then
! found element
          bonds_getind=j
          exit
         endif
        enddo
       endif
!
       end function bonds_getind
      end module bondpar
