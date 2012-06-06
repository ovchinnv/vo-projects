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
module state
 use output
!
 type varray
   real*8, pointer, dimension(:,:,:) :: v
 end type varray
!
 integer, parameter :: nvar=4
 character(len=6), parameter :: var_names(nvar) = &
& [character(len=6) :: 'EPS','KAPPA','CHARGE', 'PHI'];
 logical :: var_init(nvar)
!
 real*8, pointer, save, dimension(:,:,:) :: &
& p, & ! electrostatic potential
& eps, & ! dielectric constant
& kappa, & ! ionic strength parameter
& rhs ! charge density (source terms)
!
 logical, save :: state_initialized=.false.
!
 public state_initialize
 public state_done
!
 contains
!*******************************************************************
  subroutine state_initialize()
  use gridsize, only: nx, ny, nz, size_initialized, me, communicator
  use parser
  use output
  use files
  use formats
  implicit none
!
  character(len=16) :: whoami='STATE_INITIALIZE';
!
  type (varray), dimension(nvar) :: vars
  integer :: i, j, k, l, fid, ifmt
  real*8 :: value
  logical :: qbin
  character(len=80) :: filename, fmt, file_mode, file_format
!
  if (state_initialized) call state_done()
  if (.not. size_initialized) then
   call warning(whoami, 'SIZE NOT INITIALIZED. NOTHING DONE.', -1)
   return
  else
   allocate(p(nx,ny,nz),eps(nx,ny,nz),kappa(nx,ny,nz),rhs(nx,ny,nz))
   vars(1)%v=>eps
   vars(2)%v=>kappa
   vars(3)%v=>rhs
   vars(4)%v=>p
!
! check for assignments to eps, kappa, and charge (at present, this conflicts with 'molecule.src', in that you cannot use both assignments)
! note that the boundary points are read, too; this is necessary because interpolate epsilon & kappa values are needed at cell boundaries
!
   do i=1,nvar
! check if initialization is specified
    var_init(i)=.false.
    l=len_trim(var_names(i))
    select case(getval_nocase_upper(var_names(i)(1:l)//'INIT'))
     case('CONST','CONSTANT')
      if (existtag_nocase(var_names(i)(1:l)//'CONSTANT')) then
       vars(i)%v=atof(getval_nocase(var_names(i)(1:l)//'CONSTANT'))
       var_init(i)=.true.
      elseif (existtag_nocase(var_names(i)(1:l)//'CONST')) then
       vars(i)%v=atof(getval_nocase(var_names(i)(1:l)//'CONST'))
       var_init(i)=.true.
      else
       call warning(whoami, 'Initialization constant for variable "'//var_names(i)(1:l)//'" not specified',-1)
      endif
     case('FILE')
      if (existtag_nocase(var_names(i)(1:l)//'FILE')) then
       filename=getval_nocase(var_names(i)(1:l)//'FILE')
       call adjustleft(filename)
       k=len_trim(filename)
       call message(whoami, 'Will initialize variable "'//var_names(i)(1:l)//'" from file "'//filename(1:k)//'"')
      else
       call warning(whoami, 'Initialization file for variable "'//var_names(i)(1:l)//'" not specified',-1)
      endif
!
      if (existtag_nocase(var_names(i)(1:l)//'_MODE')) then
        file_mode=getval_nocase_upper(var_names(i)(1:l)//'_MODE')
        select case(file_mode);
         case('BINARY', 'BIN'); call message(whoami,'Using binary mode'); qbin=.true.
         case('ASCII', 'TEXT'); call message(whoami,'Using ASCII mode'); qbin=.false.
         case default ; call warning(whoami,'Unknown output mode "'//file_mode(1:len_trim(file_mode))//'" requested',-1);
         qbin=.false.
        end select
      else
        qbin=.true.
      endif
!
      ifmt=chest ! default format
      if (existtag_nocase(var_names(i)(1:l)//'_FORMAT')) then
        file_format=getval_nocase_upper(var_names(i)(1:l)//'_FORMAT');
        ifmt=-999
        do j=1,num_fmt
         if (file_format.eq.format_name(j)) then
          ifmt=j
          exit
         endif
        enddo
      endif
!
      if (ifmt.gt.0) then
       call message(whoami,format_name(ifmt)(1:len_trim(format_name(ifmt)))//' format will be used');
      else
       call warning(whoami,'Format "'//file_format(1:len_trim(file_format))//'" is not recognized',-1)
      endif
!
! if (me.le.0) &
! call message(whoami, 'Reading data for variable "'//var_names(i)(1:l)//'" from file '//filename(1:k))
!
      select case(ifmt)
       case(plot3d); call plot3Dread_scalar(filename,vars(i)%v,nx,ny,nz,1,qbin)
       case(chest); call chest_read_scalar(filename,vars(i)%v,nx,ny,nz,qbin)
      end select
!
      if (.not.fatal_warning()) var_init(i)=.true.
     case('OBJECT') ! do nothing : expect that object_grid_parameters will take care of correct initialization
      call message(whoami, 'Initialization for variable "'//var_names(i)(1:l)//'" will be computed from object')
      var_init(i)=.true.
     case default
      call warning(whoami, 'Initialization option for variable "'//var_names(i)(1:l)//'" not recognized',-1)
    end select
!
   enddo ! over all variables
!
   if (all(var_init)) state_initialized=.true.
!
  endif ! size_initialized
!
  end subroutine state_initialize
!*******************************************************************
  subroutine state_done()
  implicit none
!
  if (state_initialized) then
   deallocate(p,eps,kappa,rhs)
   state_initialized=.false.
  endif
!
  end subroutine state_done
!******************************************************************
end module state
