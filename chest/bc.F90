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
      module bc ! boundary conditions
      use parser
      use gridsize
      use files
      use constants
      private
!
!
!
      public bc_initialize
      public bc_done
      public bc_apply
      public parray
!
      type parray
       real*8, pointer, dimension(:,:) :: p
      end type parray
!
!*****************************************************************
!* boundary numbering: 1 - 'left' 2 - 'right' *
!* 3 - 'bottom' 4 - 'top' *
!* 5 - 'front' 6 - 'back' *
!*****************************************************************
      integer, parameter, public :: left=1, right=2, bottom=3, top=4, front=5, back=6
      integer, parameter, public :: west=1, east=2, south=3, north=4
! corresponding boundary names:
      character(len=4), parameter :: bc_names(6) = (/'BCX0','BCX1','BCY0','BCY1','BCZ0','BCZ1'/);
!
      integer, public :: bc_type(6) = (/0,0,0,0,0,0/);
      real*8, private :: bc_wgt(6) = (/0.,0.,0.,0.,0.,0./); ! weights that reflect the location of the boundary relative to gridpoints
!
      integer, public, parameter :: dirichlet=1, neumann=2, periodic=3, dirichletg=4, neumanng=5 ! boundary condition codes
!
      logical, public, save :: bc_initialized=.false.
      integer, public, save :: ndim=3
!
      real*8, private, pointer, dimension(:,:) :: bcx0, bcx1, bcy0, bcy1, bcz0, bcz1
      type (parray), public, dimension(6) :: bcs
!
!*************************************************************************************
      contains
       subroutine bc_initialize()
       use output
       use constants
       implicit none
!
       character(len=13) :: whoami='BC_INITIALIZE'
       character(len=80) :: filename
       character(len=20) :: keyword
!
       integer :: i, k, ii, jj, bcsize, bcsizei, bcsizej, l
       integer :: fid, iostatus
       real*8 :: value
!
       if (.not.size_initialized) then
        call warning(whoami, 'Size not initialized. Nothing done.', -1)
        return
       endif
!
       if (bc_initialized) call bc_done() ! clean up memory if needed
!
       if (q2D) then ; ndim=2 ; else ; ndim = 3 ; endif
! note: for 2d geometry, nz=3; the ghost points are ignored, hence there is only one (inner) point
! parse BC-spec from input file
       do i=1,ndim * 2
        keyword=trim(getval_nocase_upper(bc_names(i)))
        select case(keyword)
         case('DIRICHLET'); ! Dirichlet applied at boundary
          bc_type(i)=dirichlet; bc_wgt(i)=two
         case('DIRICHLETG'); ! Dirichlet applied at ghost point
          bc_type(i)=dirichletg; bc_wgt(i)=one
         case('PERIODIC');
          bc_type(i)=periodic;
         case('NEUMANN'); bc_wgt(i)=one ! valid assuming BC data have been scaled by metric
          bc_type(i)=neumann;
         case('NEUMANNG'); bc_wgt(i)=one
          bc_type(i)=neumanng; ! Neumann applied at ghostpoint
         case default
          call warning(whoami, 'Unknown boundary condition for '//bc_names(i)//'.',-1)
         return
        end select
       enddo ! over all boundaries
!
! set the actual boundary condition values
       allocate(bcx0(ny,nz), bcx1(ny,nz));
       allocate(bcy0(nx,nz), bcy1(nx,nz));
       allocate(bcz0(nx,ny), bcz1(nx,ny));
       bcs(west)%p=>bcx0; bcs(east)%p=>bcx1;
       bcs(south)%p=>bcy0; bcs(north)%p=>bcy1;
       bcs(front)%p=>bcz0; bcs(back)%p=>bcz1;
!
       do i=1, ndim * 2
        if (bc_type(i).eq.dirichlet .or. &
& bc_type(i).eq.dirichletg.or. &
& bc_type(i).eq.neumann) &
& then
         if (existtag_nocase(bc_names(i)//'CONSTANT')) then
           bcs(i)%p=atof(getval_nocase(bc_names(i)//'CONSTANT'))
         elseif (existtag_nocase(bc_names(i)//'CONST')) then
           bcs(i)%p=atof(getval_nocase(bc_names(i)//'CONST'))
         elseif (existtag_nocase(bc_names(i)//'FILE')) then
! open and read boundary condition file
          filename=getval_nocase(bc_names(i)//'FILE')
          call adjustleft(filename)
          l=len_trim(filename)
!
          bcsizei=size(bcs(i)%p,1)-2; ! size of first dimension (inner pts only)
          bcsizej=size(bcs(i)%p,2)-2; ! size of second dimension (inner pts only)
          bcsize=bcsizei*bcsizej; ! total inner size
!
          if (me.le.0) then
           call message(whoami, 'Reading boundary conditions from file '//filename(1:l))
           if (l.gt.0) then
            fid=-1
            call files_open(fid, name_=filename(1:l), form_='FORMATTED', action_='READ')
            if (fid.lt.0) call warning(whoami, 'Cannot open file. Abort.',-1)
           else
            call warning(whoami, 'File name not specified. Abort.',-1)
           endif ! l
          endif ! me
          if (fatal_warning()) return
!
          if (me.le.0) then
           k=0
           read(fid, *, IOSTAT=iostatus) value
           do while (iostatus.eq.0)
            k=k+1
! individual indices
            jj=(k-1)/bcsizei+2 ! slowly-varying second index ; offset is 2 for inner points only
            ii=mod((k-1),bcsizei)+2 ! fast-varying first index
            if (k.le.bcsize) bcs(i)%p(ii,jj)=value ! make assignment if not out of bounds of bcs(i)
            read(fid, *, IOSTAT=iostatus) value
           enddo
           call files_close(fid)
          endif ! me
! check that the right number of entries was read
!
! write(0,*) k, bcsizei, bcsizej,bcsize ! aa
          if (k.gt.bcsize) then
           call warning(whoami, 'FILE "'//filename(1:l)//'" HAS MORE ENTRIES THAN EXPECTED.', 0)
          elseif (k.lt.bcsize) then
           call warning(whoami, 'FILE "'//filename(1:l)//'" DOES NOT CONTAIN ENOUGH ENTRIES. ABORT.', -1)
           return
          endif
!
         else
           call warning(whoami, 'No values specified at boundary '//bc_names(i)//'. Abort.', -1)
           return
         endif ! read bc values
        endif ! BC data input
       enddo ! over all boundaries
!
       bc_initialized=.true.
!
       end subroutine bc_initialize
!***********************************************************************************
       subroutine bc_done
       implicit none
       integer :: i
       if (bc_initialized) then
        bc_initialized=.false.
        deallocate(bcx0,bcx1,bcy0,bcy1,bcz0,bcz1)
        do i=1,6
         nullify(bcs(i)%p)
        enddo
       endif
!
       bc_type(1:6)=0
       end subroutine bc_done
!***********************************************************************************
       subroutine bc_apply()
       use state
       use output
       implicit none
! this is the standard routine that the module should use
! it calls the utility routine below; the reason for this design is that I want to use multigrid with this code,
! which requires multiple grids
       character(len=8) :: whoami='BC_APPLY'
!
       integer :: nxm, nym, nzm
!
       interface
        subroutine apply_bc(u,g,nnx,nny,nnz,boundary,bctype,wgt)
        real*8 :: u(nnx,nny,nnz), g(*)
        integer :: nnx,nny,nnz
        integer :: boundary
        integer :: bctype
        real*8, intent(in) :: wgt
        end subroutine apply_bc
       end interface
!
       if (.not. bc_initialized) then
        call error(whoami, 'BOUNDARY CONDITIONS NOT SET UP. NOTHING DONE.', -1)
        return
       endif
!
       nxm=nx-1; nym=ny-1; nzm=nz-1;
!
       call apply_bc(p,bcs(left)% p(2:nym,2:2:nzm),nx,ny,nz,left, bc_type(left), bc_wgt(left))
       call apply_bc(p,bcs(right)% p(2:nym,2:2:nzm),nx,ny,nz,right, bc_type(right), bc_wgt(right))
       call apply_bc(p,bcs(bottom)% p(2:nxm,2:2:nzm),nx,ny,nz,bottom,bc_type(bottom),bc_wgt(bottom))
       call apply_bc(p,bcs(top)% p(2:nxm,2:2:nzm),nx,ny,nz,top, bc_type(top), bc_wgt(top))
       if (ndim.gt.2) then
        call apply_bc(p,bcs(front)% p(2:nxm,2:2:nym),nx,ny,nz,front, bc_type(front),bc_wgt(front))
        call apply_bc(p,bcs(back)% p(2:nxm,2:2:nym),nx,ny,nz,back, bc_type(back), bc_wgt(back))
       endif
!
       end subroutine bc_apply
      end module bc
!***********************************************************************************
! utility subroutine to apply boundary conditions (cannot be inside module b/c of shape casting in multigrid_solve)
      subroutine apply_bc(u,g,nnx,nny,nnz,boundary,bctype,wgt)
      use bc
      implicit none
      real*8 :: u(nnx,nny,nnz), g(*) ! the boundary array g is assumed to have ONLY inner points, but the main array u has ghost points
      real*8, intent(in) :: wgt ! bc value is determined by extrapolation using wgt, which reflects the location of the boundary
      real*8 :: wgt_
! g is assumed size recast to 1D shape because of conflict of dimension between u(3D) and g(2D) (see below)
      integer :: nnx,nny,nnz
      integer :: boundary
      integer :: bctype
!
      integer :: ib, ie, jb, je, kb, ke
      integer :: di, dj, dk, ip, jp, kp
      integer :: i, j, k, ind, dir
!
! preset indices:
!
      ib=2; ie=nnx-1; ! limiting indices of inner points
      jb=2; je=nny-1;
      kb=2; ke=nnz-1;
!
      select case(boundary)
        case(left);
         ib=ib-1; ie=ib; di=+1; dj=0; dk=0; dir=di
        case(right);
         ie=ie+1; ib=ie; di=-1; dj=0; dk=0; dir=di
        case(bottom);
         jb=jb-1; je=jb; dj=+1; di=0; dk=0; dir=dj
        case(top);
         je=je+1; jb=je; dj=-1; di=0; dk=0; dir=dj
        case(front);
         kb=kb-1; ke=kb; dk=+1; di=0; dj=0; dir=dk
         if (ndim.eq.2) ke=0; ! for two dimensions do nothing for z-boundaries
        case(back);
         ke=ke+1; kb=ke; dk=-1; di=0; dj=0; dir=dk
         if (ndim.eq.2) ke=0
        case default;
         ke=0 ! set ivalid indices so that no loops are executed; better than early return
      end select ! which boundary
!************************************
      ind=1
      select case(bctype)
!***************************
        case(dirichlet, dirichletg); ! note that array layout is crucial here: we are assuming the standard fortran organization
         do k=kb,ke; kp=k+dk;
          do j=jb,je; jp=j+dj;
           do i=ib,ie; ip=i+di;
!
            u(i,j,k) = u(ip,jp,kp) + wgt * ( g(ind) - u(ip,jp,kp) ); ! wgt: d btw gp and 1st / d btw bdry and first
!
            ind=ind+1; enddo; enddo; enddo
!***************************
        case(neumann, neumanng);
         wgt_ = wgt*(-dir)
         do k=kb,ke; kp=k+dk;
          do j=jb,je; jp=j+dj;
           do i=ib,ie; ip=i+di;
!
            u(i,j,k) = u(ip,jp,kp) + wgt_ * g(ind); ! IMPORTANT: in this case wgt is the appropriate metric
!
            ind=ind+1; enddo; enddo; enddo
!***************************
        case(periodic);
!
         dk=(nnz-2)*dk; dj=(nny-2)*dj; di=(nnx-2)*di; ! wrap around
!
         do k=kb,ke; kp=k+dk;
          do j=jb,je; jp=j+dj;
           do i=ib,ie; ip=i+di;
!
            u(i,j,k) = u(ip,jp,kp); ! g or wgt not needed at all
!***************************
                                                   enddo; enddo; enddo
      end select ! bctype
!
      end subroutine apply_bc
!***********************************************************************************
