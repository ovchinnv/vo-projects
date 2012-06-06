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
module object
 use parser
 use output
! this module will be a wrapper around various shapes that define grid properties
! in the future, may support many objects; now this seems superfluous
use molecule, only : molecule_ndim, molecule_initialize, molecule_center, molecule_align, molecule_dimens,&
molecule_done, molecule_grid_objects, molecule_ok, molecule_surface_pointer, molecule_read_parameters
 integer, parameter :: num_objects=1
 integer, parameter :: type_undefined=0
 integer, parameter :: type_molecule=1
 character(len=20), parameter, dimension(0:num_objects) :: type_names=&
  [character(len=20) :: 'UNDEFINED','MOLECULE']
!
 logical, public :: object_initialized=.false.
 integer, public :: object_type=type_undefined
!
 public object_init
 public object_read_parameters
 private object_initialize
 private object_ok
 public object_ndim
 public object_dimens
 public object_center
 public object_grid_objects
 public object_done
 public object_surface_pointer ! function that returns a pointer to the surface array in object
!
!
 contains
!*********************************************************************************!
function OBJECT_OK( )
 implicit none
 logical :: OBJECT_OK
!
 if (.not.object_initialized) call object_init()
 if (.not.object_initialized) then
  call error("OBJECT_OK" ,'NO OBJECTS ARE DEFINED.',0);
  return
 endif
! function body

!
 select case(object_type)
  case( type_molecule ); OBJECT_OK = molecule_OK( )
  case default; call error("OBJECT_OK" ,'INTERNAL ERROR: UNKNOWN OBJECT TYPE.',-1);
 end select
!

!
end function OBJECT_OK
!
!*********************************************************************************!
function OBJECT_NDIM( )
 implicit none
 integer :: OBJECT_NDIM
!
 if (.not.object_initialized) call object_init()
 if (.not.object_initialized) then
  call error("OBJECT_NDIM" ,'NO OBJECTS ARE DEFINED.',0);
  return
 endif
! function body
 OBJECT_NDIM=-1
!
 select case(object_type)
  case( type_molecule ); OBJECT_NDIM = molecule_NDIM( )
  case default; call error("OBJECT_NDIM" ,'INTERNAL ERROR: UNKNOWN OBJECT TYPE.',-1);
 end select
!

!
end function OBJECT_NDIM
!
!*********************************************************************************!
function OBJECT_SURFACE_POINTER( )
 implicit none
 real*8, pointer, dimension(:,:,:) :: OBJECT_SURFACE_POINTER
!
 if (.not.object_initialized) call object_init()
 if (.not.object_initialized) then
  call error("OBJECT_SURFACE_POINTER" ,'NO OBJECTS ARE DEFINED.',0);
  return
 endif
! function body

!
 select case(object_type)
  case( type_molecule ); OBJECT_SURFACE_POINTER => molecule_SURFACE_POINTER( )
  case default; call error("OBJECT_SURFACE_POINTER" ,'INTERNAL ERROR: UNKNOWN OBJECT TYPE.',-1);
 end select
!

!
end function OBJECT_SURFACE_POINTER
!
!********************************************************************************!
function OBJECT_DIMENS( )
 implicit none
 real*8, pointer, dimension(:,:) :: OBJECT_DIMENS
!
 if (.not.object_initialized) call object_init()
 if (.not.object_initialized) then
  call error("OBJECT_DIMENS" ,'NO OBJECTS ARE DEFINED.',0);
  return
 endif
! function body
 ALLOCATE(OBJECT_DIMENS (object_ndim(),3)); OBJECT_DIMENS=-1
!
 select case(object_type)
  case( type_molecule ); OBJECT_DIMENS = molecule_DIMENS( )
  case default; call error("OBJECT_DIMENS" ,'INTERNAL ERROR: UNKNOWN OBJECT TYPE.',-1);
 end select
!

!
end function OBJECT_DIMENS
!
!********************************************************************************!
 subroutine OBJECT_DONE( )
 implicit none

!
 if (.not.object_initialized) call object_init()
 if (.not.object_initialized) then
  call error("OBJECT_DONE" ,'NO OBJECTS ARE DEFINED.',0);
  return
 endif
! subroutine body

!
 select case(object_type)
  case( type_molecule ); CALL molecule_DONE( )
  case default; call error("OBJECT_DONE" ,'INTERNAL ERROR: UNKNOWN OBJECT TYPE.',-1);
 end select
! subroutine body
 object_initialized=.false. ; object_type=type_undefined
!
end subroutine OBJECT_DONE
!
!********************************************************************************!
 subroutine OBJECT_CENTER( rcenter, qmass )
 implicit none
 real*8 :: rcenter(3); logical :: qmass
!
 if (.not.object_initialized) call object_init()
 if (.not.object_initialized) then
  call error("OBJECT_CENTER" ,'NO OBJECTS ARE DEFINED.',0);
  return
 endif
! subroutine body

!
 select case(object_type)
  case( type_molecule ); CALL molecule_CENTER( rcenter, qmass )
  case default; call error("OBJECT_CENTER" ,'INTERNAL ERROR: UNKNOWN OBJECT TYPE.',-1);
 end select
! subroutine body

!
end subroutine OBJECT_CENTER
!
!********************************************************************************!
 subroutine OBJECT_ALIGN( qmass )
 implicit none
 logical :: qmass
!
 if (.not.object_initialized) call object_init()
 if (.not.object_initialized) then
  call error("OBJECT_ALIGN" ,'NO OBJECTS ARE DEFINED.',0);
  return
 endif
! subroutine body

!
 select case(object_type)
  case( type_molecule ); CALL molecule_ALIGN( qmass )
  case default; call error("OBJECT_ALIGN" ,'INTERNAL ERROR: UNKNOWN OBJECT TYPE.',-1);
 end select
! subroutine body

!
end subroutine OBJECT_ALIGN
!
!********************************************************************************!
 subroutine OBJECT_GRID_OBJECTS( xcen, ycen, zcen, dxcor, dycor, dzcor, eps, kappa, rhs, nx, ny, nz )
 implicit none
 real*8, dimension(:,:,:) :: eps,kappa,rhs ; real*8, dimension(:) :: xcen,ycen,zcen,dxcor,dycor,dzcor ; integer :: nx,ny,nz
!
 if (.not.object_initialized) call object_init()
 if (.not.object_initialized) then
  call error("OBJECT_GRID_OBJECTS" ,'NO OBJECTS ARE DEFINED.',0);
  return
 endif
! subroutine body

!
 select case(object_type)
  case( type_molecule ); CALL molecule_GRID_OBJECTS( xcen, ycen, zcen, dxcor, dycor, dzcor, eps, kappa, rhs, nx, ny, nz )
  case default; call error("OBJECT_GRID_OBJECTS" ,'INTERNAL ERROR: UNKNOWN OBJECT TYPE.',-1);
 end select
! subroutine body

!
end subroutine OBJECT_GRID_OBJECTS
!
!********************************************************************************!
 subroutine OBJECT_INITIALIZE( )
 implicit none

!
 if (.not.object_initialized) call object_init()
 if (.not.object_initialized) then
  call error("OBJECT_INITIALIZE" ,'NO OBJECTS ARE DEFINED.',0);
  return
 endif
! subroutine body

!
 select case(object_type)
  case( type_molecule ); CALL molecule_INITIALIZE( )
  case default; call error("OBJECT_INITIALIZE" ,'INTERNAL ERROR: UNKNOWN OBJECT TYPE.',-1);
 end select
! subroutine body

!
end subroutine OBJECT_INITIALIZE
!
!********************************************************************************!
 subroutine OBJECT_READ_PARAMETERS( )
 implicit none

!
 if (.not.object_initialized) call object_init()
 if (.not.object_initialized) then
  call error("OBJECT_READ_PARAMETERS" ,'NO OBJECTS ARE DEFINED.',0);
  return
 endif
! subroutine body

!
 select case(object_type)
  case( type_molecule ); CALL molecule_READ_PARAMETERS( )
  case default; call error("OBJECT_READ_PARAMETERS" ,'INTERNAL ERROR: UNKNOWN OBJECT TYPE.',-1);
 end select
! subroutine body

!
end subroutine OBJECT_READ_PARAMETERS
!
!********************************************************************************!
  subroutine object_init()
  implicit none
  character(len=11), parameter :: whoami='OBJECT_INIT'
  character(len=20) :: o, ou
  logical :: found, qmass, qalign
  integer :: i
  real*8 :: r_center(3)
!
  if (object_initialized) then
   call message(whoami, 'OBJECT OF TYPE '//type_names(object_type)(1:len_trim(type_names(object_type)))//' IS ALREADY INITIALIZED.')
   return
  endif
!
  if (existtag('OBJECT')) then ; o=getval('OBJECT'); ou=o; call toupper(ou);
  elseif (existtag('object')) then ; o=getval('object'); ou=o; call toupper(ou)
  else
   call warning(whoami, 'OBJECT NOT SPECIFIED.',0)
   return
  endif
! check the type of object requested and call the appropriate initialization routine
  object_type=type_undefined
!
  found=.false.
  do i=1,num_objects
   if ( ( type_names(i) .eq. o ) .or. ( type_names(i) .eq. ou ) ) then
    object_type=i
    found=.true.
    exit
   endif
  enddo
!
  if (.not.found) then
   call warning(whoami, 'OBJECT NAME '//o(1:len_trim(o))//' NOT RECOGNIZED.',0)
   return
  endif
! call auxiliary initialization routine
  object_initialized=.true. ! hack to pass initialization checking in auxiliary routine
  call object_initialize()
! check the success of the initialization :
  object_initialized=object_ok()
!
  if (.not. object_initialized) call warning(whoami,'OBJECT INITIALIZATION FAILED.',0)
!*********** process object orientation commands
  qmass=.false.
  if (existtag_nocase('OBJECT_CENTER')) then ;
   r_center=atofv(getval_nocase('OBJECT_CENTER'),3);
   if (existtag_nocase('OBJECT_MASSW')) qmass=atol(getval_nocase('OBJECT_MASSW'))
   call object_center(r_center,qmass)
  endif
!
  if (existtag_nocase('OBJECT_ALIGN')) then
   if (atol(getval_nocase('OBJECT_ALIGN'))) then
    if (existtag_nocase('OBJECT_MASSW')) qmass=atol(getval_nocase('OBJECT_MASSW'))
    call object_align(qmass);
   endif
  endif
!
  end subroutine object_init
end module object
