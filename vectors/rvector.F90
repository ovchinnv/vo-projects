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
!
      module rvector
      implicit none
!
      type real_vector
       float, dimension(:), pointer :: r
       int :: length ! length of the vector
       int :: last ! index of last element
       bool :: initialized=.false. ! has the vector been initialized
      end type real_vector
!
      private real_vector_expand
      int, parameter, private :: expand_incr=500
!
      contains
       subroutine real_vector_init( v )
       type (real_vector) :: v
! if (associated(v%r)) deallocate(v%r) ! testing unassigned pointer is an error!
       allocate(v%r(expand_incr))
       v%r=0
       v%length=expand_incr
       v%last=0
       v%initialized=.true.
       end subroutine real_vector_init
!ccccc
       subroutine real_vector_done( v )
       type (real_vector) :: v
       if (associated(v%r)) deallocate(v%r)
       v%length=0
       v%last=0
       v%initialized=.false.
       end subroutine real_vector_done
!ccccc
       subroutine real_vector_expand( v )
       type (real_vector) :: v
       int :: newlength
       float, dimension(:), allocatable :: p
!
       if (.not.v%initialized) then
        call real_vector_init(v)
       else
! assume length is valid
        newlength=v%length+expand_incr
        allocate(p(newlength)) ! new memory
        p(1:v%length)=v%r ! copy old data
        deallocate(v%r) ! delete old data
        allocate(v%r(newlength))
        v%r = p ! copy data
        deallocate(p)
        v%length=newlength
       endif
       end subroutine real_vector_expand
!ccccc
       function real_vector_add( v,i ) ! add a new element to the list (not necessarily unique)
! and return its index
       type (real_vector) :: v
       float :: i
       int :: j, real_vector_add
!
       if (.not.v%initialized) call real_vector_init(v)
! add element to the list
       if (v%last.eq.v%length) call real_vector_expand(v)
       j=v%last+1
       v%r(j)=i
       v%last=j
       real_vector_add=j
       end function real_vector_add
!ccccc
       function real_vector_uadd( v,i ) ! add a UNIQUE new element to the list and return its index
! if the element already exists, return its index
       type (real_vector) :: v
       float :: i
       int :: j, real_vector_uadd
!
       if (.not.v%initialized) call real_vector_init(v)
       do j=1,v%last
        if (v%r(j).eq.i) then
         real_vector_uadd=j
         return
        endif
       enddo
! add element to the list
       if (v%last.eq.v%length) call real_vector_expand(v)
       j=v%last+1
       v%r(j)=i
       v%last=j
       real_vector_uadd=j
       end function real_vector_uadd
!ccccc
       function real_vector_get( v,j ) ! returns v%r(j) if j is valid
       type (real_vector) :: v
       float :: i, real_vector_get
       int :: j
       if (.not.v%initialized) then
        i=-1
       elseif (j.gt.v%last.or.j.le.0) then
        i=-1
       else
        i=v%r(j)
       endif
       real_vector_get=i
       end function real_vector_get
!ccccc
       function real_vector_getlast(v) ! returns the last element, if list nonempty
       type (real_vector) :: v
       float :: i, real_vector_getlast
       if (.not.v%initialized) then
        i=-1
       elseif (v%last.le.0) then
        i=-1
       else
        i=v%r(v%last)
       endif
       real_vector_getlast=i
       end function real_vector_getlast
!ccccc
       function real_vector_getind( v,i ) result(j) ! returns j for the first "v%__DATANAME(j)=i" match
       type (real_vector) :: v
       int :: i, j
       j=-1
       if (v%initialized) then
        do j=1,v%last
         if (v%r(j).eq.i) exit
        enddo
        if (j.eq.v%last) then ; if (v%r(j).ne.i) j=-1 ; endif ! not found entry, just ran to end of loop !
       endif
       end function real_vector_getind
!ccccc
       function real_vector_delete( v,i )
       type (real_vector) :: v
       bool :: real_vector_delete
       int :: i
       if (i.gt.0.and.i.le.v%last) then ! delete
        if (i.lt.v%last) v%r(i)=v%r(v%last)
        v%last=v%last-1
        real_vector_delete=.true.
       else ! out of bounds
        real_vector_delete=.false.
       endif
       end function real_vector_delete
!ccccc
      end module rvector
