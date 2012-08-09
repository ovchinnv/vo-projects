/*#define __WRN(__WHO,__MSG) write(0,*) 'WARNING FROM: ',__WHO,': ',__MSG*/
/*#define __PRINT(__MSG) write(0,'(A)') __MSG*/
/*#define __PRINT(__MSG) call plainmessage(__MSG)*/
/*#define __PRINTL(__MSG,__LEVEL) call plainmessage(__MSG,__LEVEL)*/
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
      module atompar
! define a derived type to store atom parameters
      type atoms
       character(len=8), dimension(:), pointer :: a1
       real*8, dimension(:), pointer :: rmino2, eps, rmino214, eps14, mass
       character(len=2), dimension(:), pointer :: el ! element type (e.g. C, Na, etc)
       integer, dimension(:), pointer :: typeid ! this is the numerical ID number that appears in the topology file (and CHARMM psf)
       integer :: length ! length of the vector
       integer :: last ! index of last element
       logical :: initialized=.false. ! has the vector been initialized
      end type atoms
!
      type atom ! scalar type
       character(len=8) :: a1
       real*8 :: rmino2, eps, rmino214, eps14, mass
       character(len=2) :: el
       integer :: typeid
      end type atom
!
      integer, parameter, private :: expand_incr=200
!
      contains
!******************************************** atom data routines*********************
       subroutine atoms_init( v )
       implicit none
       type (atoms) :: v
! if (associated(v%...)) deallocate(v%...) ! testing unassigned pointer is an error!
       allocate(v%a1(expand_incr),v%rmino2(expand_incr),v%eps(expand_incr), &
& v%rmino214(expand_incr), v%eps14(expand_incr), v%mass(expand_incr),&
& v%el(expand_incr), v%typeid(expand_incr) )
!
       v%a1=''; v%rmino2=-1d0; v%rmino214=-1d0; v%eps=-1d0; v%eps14=-1d0
       v%mass=-1d0; v%el=''; v%typeid=-1
       v%length=expand_incr
       v%last=0
       v%initialized=.true.
       end subroutine atoms_init
!ccccc
       subroutine atoms_done( v )
       implicit none
       type (atoms) :: v
       if (associated(v%a1)) deallocate(v%a1)
       if (associated(v%rmino2)) deallocate(v%rmino2)
       if (associated(v%eps)) deallocate(v%eps)
       if (associated(v%rmino214)) deallocate(v%rmino214)
       if (associated(v%eps14)) deallocate(v%eps14)
       if (associated(v%mass)) deallocate(v%mass)
       if (associated(v%el)) deallocate(v%el)
       if (associated(v%typeid)) deallocate(v%typeid)
       v%length=0
       v%last=0
       v%initialized=.false.
       end subroutine atoms_done
!
       subroutine atoms_expand( v )
       implicit none
       type (atoms) :: v
       integer :: newlength
       character(len=8), dimension(:), pointer :: a1
       character(len=2), dimension(:), pointer :: el
       real*8, dimension(:), pointer :: rmino2, eps, rmino214, eps14, mass
       integer, dimension(:), pointer :: typeid
!
       if (.not.v%initialized) then
        call atoms_init(v)
       else
! assume length is valid
        newlength=v%length+expand_incr ! temporary storage space
        allocate(a1(newlength),rmino2(newlength),eps(newlength),rmino214(newlength), eps14(newlength),&
& el(newlength), mass(newlength), typeid(newlength))
        a1(1:v%length)=v%a1; rmino2(1:v%length)=v%rmino2; eps(1:v%length)=v%eps; rmino214(1:v%length)=v%rmino214
        eps14(1:v%length)=v%eps14; el(1:v%length)=v%el; mass(1:v%length)=v%mass; typeid(1:v%length)=v%typeid
        deallocate(v%a1,v%rmino2,v%eps,v%rmino214,v%eps14,v%mass,v%typeid,v%el)
        v%a1=>a1; v%rmino2=>rmino2; v%eps=>eps; v%rmino214=>rmino214; v%eps14=>eps14; v%mass=>mass; v%el=>el; v%typeid=>typeid
        v%length=newlength
       endif
       end subroutine atoms_expand
!ccccc
       function atoms_add_mass(v,a1,mass,el,typeid)
       implicit none
       type (atoms) :: v
       integer :: atoms_add_mass
       real*8 :: mass
       character(len=*) :: el
       integer :: typeid
       character(len=*) :: a1
       integer :: j
!
       if (.not.v%initialized) call atoms_init(v)
! add element to the list
       if (v%last.eq.v%length) call atoms_expand(v)
       j=v%last+1
       v%a1(j)=a1; v%mass(j)=mass; v%el(j)=el; v%typeid(j)=typeid
       v%last=j
       atoms_add_mass=j
       end function atoms_add_mass
!ccccc
       function atoms_add_nonb(v,a1,rmino2,eps,rmino214,eps14)
       implicit none
       type (atoms) :: v
       integer :: atoms_add_nonb
       real*8 :: rmino2, eps, rmino214, eps14
       integer :: typeid
       character(len=*) :: a1
       integer :: j
!
       if (.not.v%initialized) call atoms_init(v)
! add element to the list
       if (v%last.eq.v%length) call atoms_expand(v)
       j=v%last+1
       v%a1(j)=a1; v%rmino2(j)=rmino2; v%eps(j)=eps; v%rmino214(j)=rmino214; v%eps14(j)=eps14;
       v%last=j
       atoms_add_nonb=j
       end function atoms_add_nonb
!ccccc
       function atoms_uadd_mass(v,a1,mass,el,typeid)
       use output, only: warning
       implicit none
       type (atoms) :: v
       integer :: atoms_uadd_mass
       real*8 :: mass
       character(len=*) :: el
       integer :: typeid
       character(len=*) :: a1
       integer :: j
!
       if (.not.v%initialized) call atoms_init(v)
       do j=1,v%last
        if ( v%a1(j).eq.a1 ) then
! found element
         atoms_uadd_mass=j
         if (v%mass(j).gt.0.) & ! mass has been set before
          call warning('ATOMS_UADD_MASS','Parameters for atom type "'//v%a1(j)//'" already present. Will overwrite.',0)
         v%mass(j)=mass; v%el(j)=el; v%typeid(j)=typeid
         return
        endif
       enddo
! add element to the list: use regular routine
       atoms_uadd_mass=atoms_add_mass(v,a1,mass,el,typeid)
!
       end function atoms_uadd_mass
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       function atoms_uadd_nonb(v,a1,rmino2,eps,rmino214,eps14)
       use output, only: warning
       implicit none
       type (atoms) :: v
       integer :: atoms_uadd_nonb
       real*8 :: rmino2, eps, rmino214, eps14
       character(len=*) :: a1
       integer :: j
!
       if (.not.v%initialized) call atoms_init(v)
       do j=1,v%last
        if ( v%a1(j).eq.a1 ) then
! found element
         atoms_uadd_nonb=j
         if (v%rmino2(j).gt.0.) & ! vdw parms have been set before
          call warning('ATOMS_UADD_NONB','Parameters for atom type"'//v%a1(j)//'" already present. Will overwrite.',0)
         v%rmino2(j)=rmino2; v%eps(j)=eps; v%rmino214(j)=rmino214; v%eps14(j)=eps14;
         return
        endif
       enddo
! add element to the list: use regular routine
       atoms_uadd_nonb=atoms_add_nonb(v,a1,rmino2,eps,rmino214,eps14)
!
       end function atoms_uadd_nonb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       function atoms_getpar( v,a1 )
       implicit none
       type (atoms) :: v
       type (atom) :: a, atoms_getpar
       character(len=*) :: a1
       integer :: j
!
       atoms_getpar%a1='';
       atoms_getpar%typeid=-1;
!
       if (.not.v%initialized) then
!
       else
        do j=1,v%last
         if ( v%a1(j).eq.a1 ) then
! found element
          a%a1=v%a1(j); a%el=v%el(j); a%typeid=v%typeid(j);
          a%rmino2=v%rmino2(j); a%eps=v%eps(j); a%rmino214=v%rmino214(j); a%eps14=v%eps14(j); a%mass=v%mass(j)
          exit
         endif
        enddo
       endif
!
       atoms_getpar=a
       end function atoms_getpar
      end module atompar
