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
      module psfatom
! define a derived aname to store atom parameters (for PSF)
      type atomlist
       character(len=8), dimension(:), pointer :: segid, resid, resname, aname
       real*8, dimension(:), pointer :: charge, mass
       integer, dimension(:), pointer :: atomid ! psf ID (unique)
       integer, dimension(:), pointer :: typeid ! this is the numerical ID number that appears in the topology file (and CHARMM psf)
       character(len=8), dimension(:), pointer :: type ! this is the atom type (same as appears in the topology file);
       integer :: length=0 ! length of the vector
       integer :: last=0 ! index of last element
       logical :: initialized=.false. ! has the vector been initialized
      end type atomlist
!
      integer, parameter, private :: expand_incr=200
      character(len=8), parameter :: unknown='unknown?'
!
      contains
!******************************************** atom data routines *********************
       subroutine atomlist_init( v )
       implicit none
       type (atomlist) :: v
! if (associated(v%...)) deallocate(v%...) ! testing unassigned pointer is an error!
       allocate(v%segid(expand_incr),v%resid(expand_incr),v%resname(expand_incr), &
& v%aname(expand_incr), v%charge(expand_incr), v%mass(expand_incr),&
& v%typeid(expand_incr),v%type(expand_incr), v%atomid(expand_incr) )
!
       v%segid=''; v%resid=''; v%resname=''; v%aname=''; v%charge=0d0; v%type=unknown;
       v%mass=-1d0; v%typeid=-1; v%atomid=-1
       v%length=expand_incr
       v%last=0
       v%initialized=.true.
       end subroutine atomlist_init
!ccccc
       subroutine atomlist_done( v )
       implicit none
       type (atomlist) :: v
       if (associated(v%segid)) deallocate(v%segid)
       if (associated(v%resid)) deallocate(v%resid)
       if (associated(v%resname)) deallocate(v%resname)
       if (associated(v%aname)) deallocate(v%aname)
       if (associated(v%charge)) deallocate(v%charge)
       if (associated(v%mass)) deallocate(v%mass)
       if (associated(v%typeid)) deallocate(v%typeid)
       if (associated(v%type)) deallocate(v%type)
       if (associated(v%atomid)) deallocate(v%atomid)
       v%length=0
       v%last=0
       v%initialized=.false.
       end subroutine atomlist_done
!
       subroutine atomlist_expand( v )
       implicit none
       type (atomlist) :: v
       integer :: newlength
       character(len=8), dimension(:), pointer :: segid, resid, aname, resname, type
       real*8, dimension(:), pointer :: charge, mass
       integer, dimension(:), pointer :: typeid, atomid
!
       if (.not.v%initialized) then
        call atomlist_init(v)
       else
! assume length is valid
        newlength=v%length+expand_incr ! temporary storage space
        allocate(segid(newlength),resid(newlength),aname(newlength),resname(newlength), charge(newlength),&
& typeid(newlength), mass(newlength), atomid(newlength), type(newlength))
        type=unknown
        segid(1:v%length)=v%segid; resid(1:v%length)=v%resid; aname(1:v%length)=v%aname; resname(1:v%length)=v%resname
        charge(1:v%length)=v%charge; typeid(1:v%length)=v%typeid; mass(1:v%length)=v%mass; atomid(1:v%length)=v%atomid
        type(1:v%length)=v%type
        deallocate(v%segid,v%resid,v%aname,v%resname,v%charge,v%mass,v%atomid,v%typeid,v%type)
        v%segid=>segid; v%resid=>resid; v%aname=>aname; v%resname=>resname; v%charge=>charge; v%mass=>mass; v%typeid=>typeid;
        v%atomid=>atomid; v%type=>type
        v%length=newlength
       endif
       end subroutine atomlist_expand
!ccccc
       function atomlist_add(v, atomid, segid, resid, resname, aname, type, typeid, charge, mass)
! note that I am ignoring the last three entries from PSF: imove, ECH, EHA
       implicit none
       type (atomlist) :: v
       integer :: atomlist_add
       real*8 :: mass, charge
       integer :: atomid, typeid
       character(len=*) :: segid, resid, resname, aname, type
       integer :: j
!
       if (.not.v%initialized) call atomlist_init(v)
! add element to the list
       if (v%last.eq.v%length) call atomlist_expand(v)
       j=v%last+1
       v%segid(j)=segid; v%resid(j)=resid; v%resname(j)=resname;
       v%aname(j)=aname; v%mass(j)=mass; v%charge(j)=charge;
       v%typeid(j)=typeid; v%type(j)=type; v%atomid(j)=atomid
       v%last=j
       atomlist_add=j
       end function atomlist_add
!ccccc
       function atomlist_uadd_ch(v, atomid, segid, resid, resname, aname, typeid, charge, mass)
       use output, only: fout, error
       implicit none
       type (atomlist) :: v
       integer :: atomlist_uadd_ch
       real*8 :: mass, charge
       integer :: atomid, typeid
       character(len=8), parameter :: type=unknown
       character(len=*) :: segid, resid, resname, aname
       integer :: j
!
       if (.not.v%initialized) call atomlist_init(v)
       do j=1,v%last
        if ( v%atomid(j).eq.atomid ) then
! found atom
         write(fout,'(2A,I5,A)') ' ERROR(ATOMLIST_UADD): ','ATOM INDEX ',v%atomid(j),' ALREADY PRESENT. ABORT.'
         call error('ATOMLIST_UADD','',-1)
         return
        endif
       enddo
! add element to the list: use regular routine
       atomlist_uadd_ch=atomlist_add(v, atomid, segid, resid, resname, aname, type, typeid, charge, mass)
!
       end function atomlist_uadd_ch
!ccccc
       function atomlist_uadd_xp(v, atomid, segid, resid, resname, aname, type, charge, mass)
       use output, only: fout, error
       implicit none
       type (atomlist) :: v
       integer :: atomlist_uadd_xp
       real*8 :: mass, charge
       integer :: atomid
       integer, parameter :: typeid=-1
       character(len=*) :: segid, resid, resname, aname, type
       integer :: j
!
       if (.not.v%initialized) call atomlist_init(v)
       do j=1,v%last
        if ( v%atomid(j).eq.atomid ) then
! found atom
         write(fout,'(2A,I5,A)') ' ERROR(ATOMLIST_UADD): ','ATOM INDEX ',v%atomid(j),' ALREADY PRESENT. ABORT.'
         call error('ATOMLIST_UADD','',-1)
         return
        endif
       enddo
! add element to the list: use regular routine
       atomlist_uadd_xp=atomlist_add(v, atomid, segid, resid, resname, aname, type, typeid, charge, mass)
!
       end function atomlist_uadd_xp
!
      end module psfatom
