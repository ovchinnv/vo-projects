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
      module dihepar
! define a derived type to store dihe parameters
      type dihes
       character(len=8), dimension(:), pointer :: a1, a2, a3, a4
       real*8, dimension(:), pointer :: kchi, delta
       integer, dimension(:), pointer :: mult
!V(dihedral) = Kchi(1 + cos(n(chi) - delta))
!
!Kchi: kcal/mole
!n: multiplicity
!delta: degrees
       integer :: length ! length of the vector
       integer :: last ! index of last element
       logical :: initialized=.false. ! has the vector been initialized
      end type dihes
!
      type dihe
       character(len=8) :: a1, a2, a3, a4
       real*8 :: kchi, delta
       integer :: mult
      end type dihe
!
      integer, parameter, private :: expand_incr=200
!
      contains
!******************************************** dihe data routines*********************
       subroutine dihes_init( v )
       implicit none
       type (dihes) :: v
! if (associated(v%...)) deallocate(v%...) ! testing unassigned pointer is an error!
       allocate(v%a1(expand_incr),v%a2(expand_incr),v%a3(expand_incr),v%a4(expand_incr),&
                v%kchi(expand_incr),v%delta(expand_incr),v%mult(expand_incr))
!
       v%a1=''; v%a2=''; v%a3=''; v%a4=''; v%kchi=-1.; v%mult=-1; v%delta=-1.;
       v%length=expand_incr
       v%last=0
       v%initialized=.true.
       end subroutine dihes_init
!ccccc
       subroutine dihes_done( v )
       implicit none
       type (dihes) :: v
       if (associated(v%a1)) deallocate(v%a1)
       if (associated(v%a2)) deallocate(v%a2)
       if (associated(v%a3)) deallocate(v%a3)
       if (associated(v%a4)) deallocate(v%a4)
       if (associated(v%kchi)) deallocate(v%kchi)
       if (associated(v%delta)) deallocate(v%delta)
       if (associated(v%mult)) deallocate(v%mult)
       v%length=0
       v%last=0
       v%initialized=.false.
       end subroutine dihes_done
!
       subroutine dihes_expand( v )
       implicit none
       type (dihes) :: v
       integer :: newlength
       character(len=8), dimension(:), pointer :: a1, a2, a3, a4
       real*8, dimension(:), pointer :: kchi, delta
       integer, dimension(:), pointer :: mult
!
       if (.not.v%initialized) then
        call dihes_init(v)
       else
! assume length is valid
        newlength=v%length+expand_incr ! temporary storage space
        allocate(a1(newlength),a2(newlength),a3(newlength),a4(newlength),kchi(newlength),delta(newlength),mult(newlength))
        a1(1:v%length)=v%a1; a2(1:v%length)=v%a2; a3(1:v%length)=v%a3; a4(1:v%length)=v%a4; kchi(1:v%length)=v%kchi;
        mult(1:v%length)=v%mult; delta(1:v%length)=v%delta
        deallocate(v%a1,v%a2,v%a3,v%a4,v%kchi,v%delta,v%mult)
        v%a1=>a1; v%a2=>a2; v%a3=>a3; v%a4=>a4; v%kchi=>kchi; v%delta=>delta; v%mult=>mult
        v%length=newlength
       endif
       end subroutine dihes_expand
!ccccc
       function dihes_add(v,a1,a2,a3,a4,kchi,delta,mult)
       implicit none
       type (dihes) :: v
       integer :: dihes_add
       real*8 :: kchi, delta
       integer :: mult
       character(len=*) :: a1, a2, a3, a4
       integer :: j
!
       if (.not.v%initialized) call dihes_init(v)
! add element to the list
       if (v%last.eq.v%length) call dihes_expand(v)
       j=v%last+1
       v%a1(j)=a1; v%kchi(j)=kchi
       v%a2(j)=a2; v%delta(j)=delta
       v%a3(j)=a3; v%mult(j)=mult
       v%a4(j)=a4;
       v%last=j
       dihes_add=j
       end function dihes_add
!ccccc
       function dihes_uadd(v,a1,a2,a3,a4,kchi,delta,mult)
       use output, only: warning
       implicit none
       type (dihes) :: v
       integer :: dihes_uadd
       real*8 :: kchi, delta
       integer :: mult
       character(len=*) :: a1, a2, a3, a4
       integer :: j
!
       if (.not.v%initialized) call dihes_init(v)
       do j=1,v%last
        if ( v%mult(j).eq.mult.and. & ! can have several dihedral restraints with diff. multimplicities
             ((v%a1(j).eq.a1.and.v%a2(j).eq.a2.and.v%a3(j).eq.a3.and.v%a4(j).eq.a4).or.&
             (v%a1(j).eq.a4.and.v%a2(j).eq.a3.and.v%a3(j).eq.a2.and.v%a4(j).eq.a1)) ) then
! found element
         dihes_uadd=j
         call warning('DIHES_UADD','Parameters for dihe "'//v%a1(j)//'--'//v%a2(j)//'--'//v%a3(j)//'--'//v%a4(j)//&
                      '" already present. Will overwrite.',0)
         v%kchi(j)=kchi; v%delta(j)=delta; v%mult(j)=mult
! write(0,*) v%kchi(j), kchi, v%delta(j), delta, v%mult(j), mult
         return
        endif
       enddo
! add element to the list: use regular routine
       dihes_uadd=dihes_add(v,a1,a2,a3,a4,kchi,delta,mult)
!
       end function dihes_uadd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       function dihes_getind( v,a1,a2,a3,a4 )
       use ivector
       implicit none
       type (dihes) :: v
       type (int_vector) :: d
       integer, pointer :: dihes_getind(:)
       character(len=*) :: a1, a2, a3, a4
       integer :: j, k, l, m, i1, i2, flag(4)
!
       call int_vector_init(d)
       if (.not.v%initialized) then
!
       else
        do j=1,v%last
         if ( ( (v%a1(j).eq.a1.or.v%a1(j).eq.'X') &
          .and. (v%a2(j).eq.a2.or.v%a2(j).eq.'X') &
          .and. (v%a3(j).eq.a3.or.v%a3(j).eq.'X') &
          .and. (v%a4(j).eq.a4.or.v%a4(j).eq.'X') ) .or.&
              ( (v%a1(j).eq.a4.or.v%a1(j).eq.'X') &
          .and. (v%a2(j).eq.a3.or.v%a2(j).eq.'X') &
          .and. (v%a3(j).eq.a2.or.v%a3(j).eq.'X') &
          .and. (v%a4(j).eq.a1.or.v%a4(j).eq.'X') ) ) then
!
! found element
            k=int_vector_uadd(d,j)
         endif
        enddo
       endif
! now throw away redundancies (if any)
       do j=1,d%last-1
        k=j+1
        do while (k.le.d%last)
         i1=d%i(j); i2=d%i(k)
         where( (/ v%a1(i1), v%a2(i1), v%a3(i1), v%a4(i1) /).eq.'X'); flag=1; elsewhere; flag=0; endwhere; l=sum(flag);
         where( (/ v%a1(i2), v%a2(i2), v%a3(i2), v%a4(i2) /).eq.'X'); flag=1; elsewhere; flag=0; endwhere; m=sum(flag);
         if (l.gt.m) then
          d%i(j)=i2;
          d%i(k)=i1;
         else
          k=k+1;
         endif
        enddo
       enddo
!
! dihedral type must be specified using the same atom types ( am guessing that reversing atom order does not create a new type)
       do j=2, d%last
        i1=d%i(1)
        i2=d%i(j)
          if ( any( (/ v%a1(i1), v%a2(i1), v%a3(i1), v%a4(i1) /) .ne. (/ v%a1(i2), v%a2(i2), v%a3(i2), v%a4(i2) /) ).and. &
               any( (/ v%a1(i1), v%a2(i1), v%a3(i1), v%a4(i1) /) .ne. (/ v%a4(i2), v%a3(i2), v%a2(i2), v%a1(i2) /) ) &
          ) exit
       enddo
!
       allocate(dihes_getind(j))
       dihes_getind=(/j-1, d%i(1:j-1)/)
       call int_vector_done(d)
!
       end function dihes_getind
      end module dihepar
