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
      module anglpar
! define a derived type to store angle parameters
      type angles
       character(len=8), dimension(:), pointer :: a1, a2, a3
       real*8, dimension(:), pointer :: ktheta, theta0, kub, s0
! V(angle) = Kb(b - b0)**2
!Kb: kcal/mole/A**2
!b0: A
       integer :: length ! length of the vector
       integer :: last ! index of last element
       logical :: initialized=.false. ! has the vector been initialized
      end type angles
!
      integer, parameter, private :: expand_incr=200
!
      contains
!******************************************** angle data routines*********************
       subroutine angles_init( v )
       implicit none
       type (angles) :: v
! if (associated(v%...)) deallocate(v%...) ! testing unassigned pointer is an error!
       allocate(v%a1(expand_incr),v%a2(expand_incr),v%a3(expand_incr),v%ktheta(expand_incr),v%theta0(expand_incr),&
& v%kub(expand_incr),v%s0(expand_incr))
!
       v%a1=''; v%a2=''; v%a3=''; v%ktheta=0; v%theta0=0; v%kub=0; v%s0=0
       v%length=expand_incr
       v%last=0
       v%initialized=.true.
       end subroutine angles_init
!ccccc
       subroutine angles_done( v )
       implicit none
       type (angles) :: v
       if (associated(v%a1)) deallocate(v%a1)
       if (associated(v%a2)) deallocate(v%a2)
       if (associated(v%a3)) deallocate(v%a3)
       if (associated(v%ktheta)) deallocate(v%ktheta)
       if (associated(v%theta0)) deallocate(v%theta0)
       if (associated(v%kub)) deallocate(v%kub)
       if (associated(v%s0)) deallocate(v%s0)
       v%length=0
       v%last=0
       v%initialized=.false.
       end subroutine angles_done
!
       subroutine angles_expand( v )
       implicit none
       type (angles) :: v
       integer :: newlength
       character(len=8), dimension(:), pointer :: a1, a2, a3
       real*8, dimension(:), pointer :: kub, s0, ktheta, theta0
!
       if (.not.v%initialized) then
        call angles_init(v)
       else
! assume length is valid
        newlength=v%length+expand_incr ! temporary storage space
        allocate(a1(newlength),a2(newlength),a3(newlength),kub(newlength),s0(newlength),&
& ktheta(newlength),theta0(newlength)) ! copy old data
        a1(1:v%length)=v%a1; a2(1:v%length)=v%a2; a3(1:v%length)=v%a3; kub(1:v%length)=v%kub; s0(1:v%length)=v%s0
        ktheta(1:v%length)=v%ktheta; theta0(1:v%length)=v%theta0
        deallocate(v%a1,v%a2,v%a3,v%ktheta,v%theta0,v%kub,v%s0)
        v%a1=>a1; v%a2=>a2; v%a3=>a3; v%ktheta=>ktheta; v%theta0=>theta0; v%kub=>kub; v%s0=>s0;
        v%length=newlength
       endif
       end subroutine angles_expand
!ccccc
       function angles_add(v,a1,a2,a3,ktheta,theta0,kub,s0)
       implicit none
       type (angles) :: v
       integer :: angles_add
       real*8 :: ktheta, theta0, kub, s0
       character(len=*) :: a1, a2, a3
       integer :: j
!
       if (.not.v%initialized) call angles_init(v)
! add element to the list
       if (v%last.eq.v%length) call angles_expand(v)
       j=v%last+1
       v%a1(j)=a1; v%a2(j)=a2; v%a3(j)=a3; v%ktheta(j)=ktheta; v%theta0(j)=theta0; v%kub(j)=kub; v%s0(j)=s0
       v%last=j
       angles_add=j
       end function angles_add
!ccccc
       function angles_uadd(v,a1,a2,a3,ktheta,theta0,kub,s0)
       use output, only: warning
       implicit none
       type (angles) :: v
       integer :: angles_uadd
       real*8 :: ktheta, theta0, kub, s0
       character(len=*) :: a1, a2, a3
       integer :: j
!
       if (.not.v%initialized) call angles_init(v)
       do j=1,v%last
        if ( v%a2(j).eq.a2.and.((v%a1(j).eq.a1.and.v%a3(j).eq.a3).or.(v%a1(j).eq.a3.and.v%a3(j).eq.a1)) ) then
! found element
         angles_uadd=j
         call warning('ANGLES_UADD','Parameters for angle "'//a1//'--'//a2//'--'//a3//'" already present. Will overwrite.',0)
         v%ktheta(j)=ktheta; v%theta0(j)=theta0; v%kub(j)=kub; v%s0(j)=s0
         return
        endif
       enddo
! add element to the list: use regular routine
       angles_uadd=angles_add(v,a1,a2,a3,ktheta,theta0,kub,s0)
!
       end function angles_uadd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       function angles_getind( v,a1,a2,a3 )
       implicit none
       type (angles) :: v
       real*8 :: angles_getind
       character(len=*) :: a1, a2, a3
       integer :: j
!
       angles_getind=-999
       if (.not.v%initialized) then
!
       else
        do j=1,v%last
         if ( v%a2(j).eq.a2.and.((v%a1(j).eq.a1.and.v%a3(j).eq.a3).or.(v%a1(j).eq.a3.and.v%a3(j).eq.a1)) ) then
! found element
          angles_getind=j
          exit
         endif
        enddo
       endif
!
       end function angles_getind
      end module anglpar
