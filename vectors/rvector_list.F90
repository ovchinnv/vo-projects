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
      module rvector_list
      use rvector
      implicit none
!
      type real_vlist
       int, dimension(:), pointer :: i ! int label
       type (real_vector), dimension(:), pointer :: v ! list of vectors
       int :: length ! length of the vector
       int :: last ! index of last element
       bool :: initialized=.false. ! has the vector been initialized
      end type real_vlist
!
      private real_vlist_expand
      private real_vlist_add
      int, parameter, private :: expand_incr=1
!
      contains
       subroutine real_vlist_init( vl, length )
       int, optional :: length
       int :: i, llength
       type (real_vlist) :: vl
! if (associated(vl%v)) then ... ! testing unassigned pointer is an error in F90 !
       if (present(length)) then ; llength=length;
       else ; llength=expand_incr ; endif
       allocate(vl%v(llength))
       allocate(vl%i(llength))
       vl%length=llength
       vl%last=0
! allocate new vectors and point to them
       do i=1,vl%length ; call real_vector_init(vl%v(i)) ; enddo
       vl%i=0
       vl%initialized=.true.
       end subroutine real_vlist_init
!
       subroutine real_vlist_done( vl )
       int :: i
       type (real_vlist) :: vl
       if (associated(vl%v)) then
        do i=1,vl%length ; call real_vector_done(vl%v(i)) ; enddo
        deallocate(vl%v)
        vl%length=0
        vl%last=0
       endif
       if (associated(vl%i)) deallocate(vl%i)
       vl%initialized=.false.
       end subroutine real_vlist_done
!
       subroutine real_vlist_expand( vl )
       type (real_vlist) :: vl
       int :: newlength, i
!
       type (real_vlist) :: wl ! temporary list
!
       if (.not.vl%initialized) then
        call real_vlist_init(vl)
       else
        newlength=vl%length+expand_incr
!
        allocate(wl%i(newlength)) ! new memory for labels
        wl%i(1:vl%length)=vl%i ! copy old data
        deallocate(vl%i) ! delete old data
        allocate(vl%i(newlength))
        vl%i = wl%i ! copy data
!
        allocate(wl%v(vl%length))
        do i=1,vl%length ; wl%v(i)=vl%v(i) ; enddo ! do not quite know why this should work in F90, but it seems to ...
        deallocate(vl%v) ! delete old data
        allocate(vl%v(newlength))
        do i=1,vl%length ; vl%v(i)=wl%v(i) ; enddo ! copy existing vectors
        do i=vl%length+1,newlength;call real_vector_init(vl%v(i));enddo ! allocate space for additional vectors
        deallocate(wl%v)
        vl%length=newlength
       endif
       end subroutine real_vlist_expand
!
       function real_vlist_add( vl, i, j ) ! add a new list labeled 'i' ; then add element 'j' to it & return index of label i
       type (real_vlist) :: vl
       int :: i, k, l, real_vlist_add
       float, optional :: j
!
       if (.not.vl%initialized) call real_vlist_init(vl)
       if (vl%last.eq.vl%length) call real_vlist_expand(vl)
       k=vl%last+1
       vl%i(k)=i ! new list label
       if (present(j)) l=real_vector_add(vl%v(k),j) ! add element to the list
       vl%last=k
       real_vlist_add=k
       end function real_vlist_add
!
       function real_vlist_uadd( vl, i, j ) ! add a new element to the list and return the index of the list
       type (real_vlist) :: vl
       int :: i, k, l, real_vlist_uadd
       float, optional :: j
!
       if (.not.vl%initialized) call real_vlist_init(vl)
       do k=1,vl%last
        if (vl%i(k).eq.i) then
         real_vlist_uadd=k
         if (present(j)) l=real_vector_add(vl%v(k),j)
         return
        endif
       enddo
! if we are here, that means a matching entry was not found; we can safely call regular add (not unique)
       if (present(j)) then
        real_vlist_uadd=real_vlist_add( vl, i, j )
       else
        real_vlist_uadd=real_vlist_add( vl, i )
       endif
!
       end function real_vlist_uadd
!
       function real_vlist_uaddu( vl, i, j ) ! add a UNIQUE new element to the list and return the index of the list
       type (real_vlist) :: vl
       int :: i, k, l, real_vlist_uaddu
       float, optional :: j
!
       if (.not.vl%initialized) call real_vlist_init(vl)
       do k=1,vl%last
        if (vl%i(k).eq.i) then
         real_vlist_uaddu=k
         if (present(j)) l=real_vector_uadd(vl%v(k),j)
         return
        endif
       enddo
! if we are here, that means a matching entry was not found; we can safely call regular add
       if (present(j)) then
        real_vlist_uaddu=real_vlist_add( vl, i, j )
       else
        real_vlist_uaddu=real_vlist_add( vl, i )
       endif
!
       end function real_vlist_uaddu
!
       function real_vlist_get( vl, i ) ! if i is found, returns the corresponding list
       type (real_vlist) :: vl
       int :: i, k
       float, pointer :: real_vlist_get(:)
       do k=1,vl%last
         if (vl%i(k).eq.i) then
          allocate(real_vlist_get(vl%v(k)%last))
          real_vlist_get=vl%v(k)%r(1:vl%v(k)%last)
          return
         endif
       enddo
       end function real_vlist_get
!
       function real_vlist_getind( vl,i ) result(j) ! returns j for the first "vl%i(j)=i" match
       type (real_vlist) :: vl
       int :: i, j
       j=-1
       if (vl%initialized) then
        do j=1,vl%last
         if (vl%i(j).eq.i) exit
        enddo
        if (j.eq.vl%last) then ; if (vl%i(j).ne.i) j=-1 ; endif ! not found entry, just ran to end of loop !
       endif
       end function real_vlist_getind
!
       function real_vlist_delete( vl,i ) ! delete list that corresponds to the tag 'i'
       type (real_vlist) :: vl
       bool :: real_vlist_delete
       int :: i, k, l
       real_vlist_delete=.false.
       do k=1,vl%last
         if (vl%i(k).eq.i) then
          l=vl%last
          vl%i(k)=vl%i(l)
          call real_vector_done(vl%v(k))
          vl%v(k)=vl%v(l)
          vl%last=vl%last-1
! nullify pointer at position l; otherwise it will hang when we deallocate position i
          nullify(vl%v(l)%r)
          real_vlist_delete=.true.
          exit
         endif
       enddo
!
       end function real_vlist_delete
!
      end module rvector_list
