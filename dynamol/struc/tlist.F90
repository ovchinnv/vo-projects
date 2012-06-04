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
      module tlist
! define a derived type to store topology lists (atom, bond, angle, dihe)
      type toplist
       integer, dimension(:,:), pointer :: ind
       integer :: length ! length of the vector
       integer :: last ! index of last element
       logical :: initialized=.false. ! has the vector been initialized
      end type toplist
!
      integer, parameter, private :: expand_incr=200
!
      contains
!******************************************** angle data routines*********************
       subroutine toplist_init( v, d )
       implicit none
       type (toplist) :: v
       integer :: d
! if (associated(v%...)) deallocate(v%...) ! testing unassigned pointer is an error!
       allocate(v%ind(d,expand_incr))
       v%ind=-1;
       v%length=expand_incr
       v%last=0
       v%initialized=.true.
       end subroutine toplist_init
!ccccc
       subroutine toplist_done( v )
       implicit none
       type (toplist) :: v
       if (associated(v%ind)) deallocate(v%ind)
       v%length=0
       v%last=0
       v%initialized=.false.
       end subroutine toplist_done
!
       subroutine toplist_expand( v,d )
       implicit none
       type (toplist) :: v
       integer :: d, newlength
       integer, dimension(:,:), pointer :: ind
!
       if (.not.v%initialized) then
        call toplist_init(v,d)
       else
! assume length is valid
        newlength=v%length+expand_incr
        allocate(ind(d,newlength))
        ind(:,1:v%length)=v%ind;
        deallocate(v%ind)
        v%ind=>ind
        v%length=newlength
       endif
       end subroutine toplist_expand
!ccccc
       function toplist_add(v,d,ind)
       implicit none
       type (toplist) :: v
       integer :: toplist_add
       integer :: d
       integer :: ind(d)
       integer :: j
!
       if (.not.v%initialized) call toplist_init(v,d)
! add element to the list
       if (v%last.eq.v%length) call toplist_expand(v,d)
       j=v%last+1
       v%ind(:,j)=ind;
       v%last=j
       toplist_add=j
       end function toplist_add
!ccccc
       function toplist_uadd(v,d,ind)
       use output, only: fout
       implicit none
       type (toplist) :: v
       integer :: toplist_uadd
       integer :: d
       integer :: ind(d)
       integer :: j
       character(len=10) :: fmt
!
       if (.not.v%initialized) call toplist_init(v,d)
       do j=1,v%last
        if (all(v%ind(1:d,j).eq.ind).or.all(v%ind(1:d,j).eq.ind(d:1:-1))) then ! assuming that reversing order generates an equivalent entry
! found element
         write(fmt,'(I10)') d
         write(fout,'(2A,'//fmt//'I5,A)') ' NONFATAL WARNING (TOPLIST_UADD)','Topology entry ',ind,'already present.' ! this can be OK -- just warn
         exit
        endif
       enddo
! add element to the list: use regular routine
       toplist_uadd=toplist_add(v,d,ind)
!
       end function toplist_uadd
      end module tlist
