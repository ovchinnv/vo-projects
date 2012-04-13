! This source file was was generated automatically from a master source
! code tree, which may or may not be distributed with this code,
! because it is up to the distributor, and not up to me.
! If you edit this file (rather than the master source file)
! your changes will be lost if another pull from the master tree occurs.
! In case you are wondering why, this approach makes it possible for
! me to have the same master source code interfaced with different
! applications (some of which are written in a way that is very far
! from being object-oriented) at the source level
!
!
      module ivector
      implicit none
!
      type int_vector
       int, dimension(:), pointer :: i
       int :: length ! length of the vector
       int :: last ! index of last element
       bool :: initialized=.false. ! has the vector been initialized
      end type int_vector
!
      private int_vector_expand
      int, parameter, private :: expand_incr=500
!
      contains
       subroutine int_vector_init( v )
       type (int_vector) :: v
! if (associated(v%i)) deallocate(v%i) ! testing unassigned pointer is an error!
       allocate(v%i(expand_incr))
       v%i=0
       v%length=expand_incr
       v%last=0
       v%initialized=.true.
       end subroutine int_vector_init
!ccccc
       subroutine int_vector_done( v )
       type (int_vector) :: v
       if (associated(v%i)) deallocate(v%i)
       v%length=0
       v%last=0
       v%initialized=.false.
       end subroutine int_vector_done
!ccccc
       subroutine int_vector_expand( v )
       type (int_vector) :: v
       int :: newlength
       int, dimension(:), allocatable :: p
!
       if (.not.v%initialized) then
        call int_vector_init(v)
       else
! assume length is valid
        newlength=v%length+expand_incr
        allocate(p(newlength)) ! new memory
        p(1:v%length)=v%i ! copy old data
        deallocate(v%i) ! delete old data
        allocate(v%i(newlength))
        v%i = p ! copy data
        deallocate(p)
        v%length=newlength
       endif
       end subroutine int_vector_expand
!ccccc
       function int_vector_add( v,i ) ! add a new element to the list (not necessarily unique)
! and return its index
       type (int_vector) :: v
       int :: i
       int :: j, int_vector_add
!
       if (.not.v%initialized) call int_vector_init(v)
! add element to the list
       if (v%last.eq.v%length) call int_vector_expand(v)
       j=v%last+1
       v%i(j)=i
       v%last=j
       int_vector_add=j
       end function int_vector_add
!ccccc
       function int_vector_uadd( v,i ) ! add a UNIQUE new element to the list and return its index
! if the element already exists, return its index
       type (int_vector) :: v
       int :: i
       int :: j, int_vector_uadd
!
       if (.not.v%initialized) call int_vector_init(v)
       do j=1,v%last
        if (v%i(j).eq.i) then
         int_vector_uadd=j
         return
        endif
       enddo
! add element to the list
       if (v%last.eq.v%length) call int_vector_expand(v)
       j=v%last+1
       v%i(j)=i
       v%last=j
       int_vector_uadd=j
       end function int_vector_uadd
!ccccc
       function int_vector_get( v,j ) ! returns v%i(j) if j is valid
       type (int_vector) :: v
       int :: i, int_vector_get
       int :: j
       if (.not.v%initialized) then
        i=-1
       elseif (j.gt.v%last.or.j.le.0) then
        i=-1
       else
        i=v%i(j)
       endif
       int_vector_get=i
       end function int_vector_get
!ccccc
       function int_vector_getlast(v) ! returns the last element, if list nonempty
       type (int_vector) :: v
       int :: i, int_vector_getlast
       if (.not.v%initialized) then
        i=-1
       elseif (v%last.le.0) then
        i=-1
       else
        i=v%i(v%last)
       endif
       int_vector_getlast=i
       end function int_vector_getlast
!ccccc
       function int_vector_getind( v,i ) result(j) ! returns j for the first "v%__DATANAME(j)=i" match
       type (int_vector) :: v
       int :: i, j
       j=-1
       if (v%initialized) then
        do j=1,v%last
         if (v%i(j).eq.i) exit
        enddo
        if (j.eq.v%last) then ; if (v%i(j).ne.i) j=-1 ; endif ! not found entry, just ran to end of loop !
       endif
       end function int_vector_getind
!ccccc
       function int_vector_delete( v,i )
       type (int_vector) :: v
       bool :: int_vector_delete
       int :: i
       if (i.gt.0.and.i.le.v%last) then ! delete
        if (i.lt.v%last) v%i(i)=v%i(v%last)
        v%last=v%last-1
        int_vector_delete=.true.
       else ! out of bounds
        int_vector_delete=.false.
       endif
       end function int_vector_delete
!ccccc
      end module ivector
