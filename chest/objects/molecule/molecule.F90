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
module molecule
 use system
 use output
 use parser
 use constants
!
 logical, private :: molecule_initialized=.false.
 real*8, private, allocatable :: epsatom(:) ! relative dielectric as a function of atom (not currently used)
!
 real*8, private :: cutoff_eps_stdev, cutoff_eps ! grid smoothing cutoff normalized by standard deviation of gaussian (epsilon)
 real*8, private :: cutoff_kappa_stdev, cutoff_kappa ! grid smoothing cutoff normalized by standard deviation of gaussian (kappa)
 real*8, private :: cutoff_charge_stdev, cutoff_charge ! grid smoothing cutoff normalized by standard deviation of gaussian (charges)
 real*8, private :: cutoff_surf_stdev, cutoff_surf ! grid smoothing cutoff normalized by standard deviation of gaussian (surface)
!
 real*8, private :: stdev_eps, oo_stdev_eps! , oo_stdev_eps2 ! std. dev. of the gaussian filter (epsilon)
 real*8, private :: stdev_kappa, oo_stdev_kappa! , oo_stdev_kappa2
 real*8, private :: stdev_charge, oo_stdev_charge! , oo_stdev_charge2
 real*8, private :: stdev_surf, oo_stdev_surf! , oo_stdev_surf2
!
 real*8, private :: eps_padding, kappa_padding, charge_padding, surf_padding ! padding values to radii which might account for e.g. solvent layer(s)
!
 real*8, private :: eps_solute, eps_solvent ! constant values
 real*8, private :: kappa_solute, kappa_solvent ! ionic strength parameter normally, kappa_solute zero
!
 integer, private, allocatable :: closest_pt(:,:)
 real*8, public, pointer :: rho(:,:,:) ! surface 'density', the isosurfaces of which are used to assign the molecular surface
!
 real*8, private, parameter :: osq2=1d0/sqrt(2d0)
 real*8, private, parameter :: osq2pi=1d0/sqrt(twopi)
!
!************************************** DEFAULT VALUES
!
 real*8, private, parameter :: eps_solute_default = 2.5d0,&
                              eps_solvent_default = 80d0
 real*8, private, parameter :: kappa_solvent_default = 0.15d0,&
                              kappa_solute_default= 0d0
 real*8, private, parameter :: cutoff_stdev_default =3d0, & ! dimensionless (see above)
                              padding_default =0d0, & ! units of radius
                              stdev_default = - 1d0 ! MEANS: determine from grid in grid_data
!
 real*8, private, parameter :: cutoff_eps_stdev_default = cutoff_stdev_default,&
                              cutoff_kappa_stdev_default = cutoff_stdev_default,&
                              cutoff_charge_stdev_default= cutoff_stdev_default,&
                              cutoff_surf_stdev_default = cutoff_stdev_default
!
 real*8, private, parameter :: stdev_eps_default = stdev_default,&
                              stdev_kappa_default = stdev_default,&
                              stdev_charge_default = stdev_default,&
                              stdev_surf_default = stdev_default
!
!
 real*8, private, parameter :: kappa_padding_default = padding_default,&
                              charge_padding_default = padding_default,&
                              eps_padding_default = padding_default,&
                              surf_padding_default = padding_default
! define character variables for output
! preprocess with a c++ - compatible preprocessor
 character(len=20), private :: eps_solute_default_str
 character(len=20), private :: eps_solvent_default_str
 character(len=20), private :: kappa_solute_default_str
 character(len=20), private :: kappa_solvent_default_str
 character(len=20), private :: stdev_eps_default_str
 character(len=20), private :: stdev_kappa_default_str
 character(len=20), private :: stdev_charge_default_str
 character(len=20), private :: stdev_surf_default_str
 character(len=20), private :: cutoff_eps_stdev_default_str
 character(len=20), private :: cutoff_kappa_stdev_default_str
 character(len=20), private :: cutoff_charge_stdev_default_str
 character(len=20), private :: cutoff_surf_stdev_default_str
 character(len=20), private :: kappa_padding_default_str
 character(len=20), private :: eps_padding_default_str
 character(len=20), private :: charge_padding_default_str
 character(len=20), private :: surf_padding_default_str
!
!***************************** subroutines ****************
 public molecule_initialize ! initialize molecular coordinates, radii and charges
 public molecule_ok
 public molecule_read_parameters ! query parser for molecule parameters
 public molecule_center ! translate center to origin
 public molecule_align ! align principal vectors with coordinate system
 public molecule_dimens ! translate to center
 public molecule_done ! deallocate memory
 public molecule_grid_objects ! public routine used to obtain grid values from this module
 public molecule_ndim ! return number of dimensions (inherited from system)
 public molecule_surface_pointer
!
!**********************************************************
 contains
!********************************************************************
  logical function molecule_ok(); molecule_ok=molecule_initialized;end function molecule_ok
!********************************************************************
  function molecule_surface_pointer()
  implicit none
  real*8, pointer, dimension(:,:,:) :: molecule_surface_pointer
  nullify(molecule_surface_pointer)
  if (molecule_initialized) then
   if (associated(rho)) molecule_surface_pointer=>rho
  endif
  end function molecule_surface_pointer
!********************************************************************
  subroutine molecule_read_parameters
  implicit none
  integer :: l
  character(len=24), parameter :: whoami='MOLECULE_READ_PARAMETERS'
  character(len=10) :: keyword
! character(len=7), parameter :: fmt='(F10.5)'
  character(len=15), parameter :: fmt='('//realfmt//')'
! convert default parameter values to character arrays
 write( eps_solute_default_str , fmt) eps_solute_default
 write( eps_solvent_default_str , fmt) eps_solvent_default
 write( kappa_solute_default_str , fmt) kappa_solute_default
 write( kappa_solvent_default_str , fmt) kappa_solvent_default
 write( cutoff_eps_stdev_default_str , fmt) cutoff_eps_stdev_default
 write( cutoff_kappa_stdev_default_str , fmt) cutoff_kappa_stdev_default
 write( cutoff_charge_stdev_default_str , fmt) cutoff_charge_stdev_default
 write( cutoff_surf_stdev_default_str , fmt) cutoff_surf_stdev_default
 write( stdev_eps_default_str , fmt) stdev_eps_default
 write( stdev_kappa_default_str , fmt) stdev_kappa_default
 write( stdev_charge_default_str , fmt) stdev_charge_default
 write( stdev_surf_default_str , fmt) stdev_surf_default
 write( eps_padding_default_str , fmt) eps_padding_default
 write( kappa_padding_default_str , fmt) kappa_padding_default
 write( charge_padding_default_str , fmt) charge_padding_default
 write( surf_padding_default_str , fmt) surf_padding_default
!
! query parser for parameters
!
! solvent dielectric
if (.not.existtag_nocase("eps_solvent")) then
   call warning(whoami, "Solvent dielectric"//' unspecified. Defaulting to '// &
                         eps_solvent_default_str// '.',0)
   eps_solvent=eps_solvent_default
  else
   keyword=getval("eps_solvent");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
                        "Solvent dielectric"//' to '//keyword(1:l)//'.');
   eps_solvent = atof(keyword)
   if (eps_solvent .lt. 0d0) &
   call warning(whoami,"Solvent dielectric"// ' < 0 (Is this what you want?)',0)
endif
!#undef
!#undef __COMPUTE
!
! solute dielectric
if (.not.existtag_nocase("eps_solute")) then
   call warning(whoami, "Solute dielectric"//' unspecified. Defaulting to '// &
                         eps_solute_default_str// '.',0)
   eps_solute=eps_solute_default
  else
   keyword=getval("eps_solute");
   select case(keyword)
    case('PDB','PQR','CHARMM','ATOMID')
    call warning(whoami, 'Atom-based parameters not yet supported. Abort.',-1)
    return
   end select
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
                        "Solute dielectric"//' to '//keyword(1:l)//'.');
   eps_solute = atof(keyword)
   if (eps_solute .lt. 0d0) &
   call warning(whoami,"Solute dielectric"// ' < 0 (Is this what you want?)',0)
endif
!#undef
!#undef __COMPUTE
!
! solvent ionic strength
if (.not.existtag_nocase("kappa_solvent")) then
   call warning(whoami, "Solvent ionic strength"//' unspecified. Defaulting to '// &
                         kappa_solvent_default_str// '.',0)
   kappa_solvent=kappa_solvent_default
  else
   keyword=getval("kappa_solvent");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
                        "Solvent ionic strength"//' to '//keyword(1:l)//'.');
   kappa_solvent = atof(keyword)
   if (kappa_solvent .lt. 0d0) &
   call warning(whoami,"Solvent ionic strength"// ' < 0 (Is this what you want?)',0)
endif
!#undef
!#undef __COMPUTE
! solute ionic strength
if (.not.existtag_nocase("kappa_solute")) then
   call warning(whoami, "Solute ionic strength"//' unspecified. Defaulting to '// &
                         kappa_solute_default_str// '.',0)
   kappa_solute=kappa_solute_default
  else
   keyword=getval("kappa_solute");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
                        "Solute ionic strength"//' to '//keyword(1:l)//'.');
   kappa_solute = atof(keyword)
   if (kappa_solute .lt. 0d0) &
   call warning(whoami,"Solute ionic strength"// ' < 0 (Is this what you want?)',0)
endif
!#undef
!#undef __COMPUTE
!%%%%%%%%%%%%%%%%%%%% padding parameters for grid smoothing %%%%%%%%%%%%%%%%%%%%%%%%%%%%
! charge padding
if (.not.existtag_nocase("charge_padding")) then
   call warning(whoami, "Charge padding"//' unspecified. Defaulting to '// &
                         charge_padding_default_str// '.',0)
   charge_padding=charge_padding_default
  else
   keyword=getval("charge_padding");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
                        "Charge padding"//' to '//keyword(1:l)//'.');
   charge_padding = atof(keyword)
   if ( charge_padding .lt. 0d0) &
   call warning(whoami,"Charge padding"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef
!#undef __COMPUTE
! kappa padding
if (.not.existtag_nocase("kappa_padding")) then
   call warning(whoami, "Ionic strength padding"//' unspecified. Defaulting to '// &
                         kappa_padding_default_str// '.',0)
   kappa_padding=kappa_padding_default
  else
   keyword=getval("kappa_padding");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
                        "Ionic strength padding"//' to '//keyword(1:l)//'.');
   kappa_padding = atof(keyword)
   if ( kappa_padding .lt. 0d0) &
   call warning(whoami,"Ionic strength padding"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef
!#undef __COMPUTE
! epsilon padding
if (.not.existtag_nocase("eps_padding")) then
   call warning(whoami, "Dielectric padding"//' unspecified. Defaulting to '// &
                         eps_padding_default_str// '.',0)
   eps_padding=eps_padding_default
  else
   keyword=getval("eps_padding");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
                        "Dielectric padding"//' to '//keyword(1:l)//'.');
   eps_padding = atof(keyword)
   if ( eps_padding .lt. 0d0) &
   call warning(whoami,"Dielectric padding"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef
!#undef __COMPUTE
! density padding
if (.not.existtag_nocase("surf_padding")) then
   call warning(whoami, "Density padding"//' unspecified. Defaulting to '// &
                         surf_padding_default_str// '.',0)
   surf_padding=surf_padding_default
  else
   keyword=getval("surf_padding");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
                        "Density padding"//' to '//keyword(1:l)//'.');
   surf_padding = atof(keyword)
   if ( surf_padding .lt. 0d0) &
   call warning(whoami,"Density padding"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef
!#undef __COMPUTE
!%%%%%%%%%%%%%% cutoffs beyond which the grid values are assumed to be zero %%%%%%%%%%%%%%%%%%%%%%
! cutoff for dielectric
if (.not.existtag_nocase("cutoff_eps_stdev")) then
   call warning(whoami, "Cutoff for dielectric smoothing"//' unspecified. Defaulting to '// &
                         cutoff_eps_stdev_default_str// '.',0)
   cutoff_eps_stdev=cutoff_eps_stdev_default
  else
   keyword=getval("cutoff_eps_stdev");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
                        "Cutoff for dielectric smoothing"//' to '//keyword(1:l)//'.');
   cutoff_eps_stdev = atof(keyword)
   if ( cutoff_eps_stdev .lt. 0d0) &
   call warning(whoami,"Cutoff for dielectric smoothing"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef
!#undef __COMPUTE
! cutoff for ionic strength
if (.not.existtag_nocase("cutoff_kappa_stdev")) then
   call warning(whoami, "Cutoff for ionic strength smoothing"//' unspecified. Defaulting to '// &
                         cutoff_kappa_stdev_default_str// '.',0)
   cutoff_kappa_stdev=cutoff_kappa_stdev_default
  else
   keyword=getval("cutoff_kappa_stdev");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
                        "Cutoff for ionic strength smoothing"//' to '//keyword(1:l)//'.');
   cutoff_kappa_stdev = atof(keyword)
   if ( cutoff_kappa_stdev .lt. 0d0) &
   call warning(whoami,"Cutoff for ionic strength smoothing"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef
!#undef __COMPUTE
! cutoff for charge
if (.not.existtag_nocase("cutoff_charge_stdev")) then
   call warning(whoami, "Cutoff for charge smoothing"//' unspecified. Defaulting to '// &
                         cutoff_charge_stdev_default_str// '.',0)
   cutoff_charge_stdev=cutoff_charge_stdev_default
  else
   keyword=getval("cutoff_charge_stdev");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
                        "Cutoff for charge smoothing"//' to '//keyword(1:l)//'.');
   cutoff_charge_stdev = atof(keyword)
   if ( cutoff_charge_stdev .lt. 0d0) &
   call warning(whoami,"Cutoff for charge smoothing"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef
!#undef __COMPUTE
! cutoff for density
if (.not.existtag_nocase("cutoff_surf_stdev")) then
   call warning(whoami, "Cutoff for density smoothing"//' unspecified. Defaulting to '// &
                         cutoff_surf_stdev_default_str// '.',0)
   cutoff_surf_stdev=cutoff_surf_stdev_default
  else
   keyword=getval("cutoff_surf_stdev");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
                        "Cutoff for density smoothing"//' to '//keyword(1:l)//'.');
   cutoff_surf_stdev = atof(keyword)
   if ( cutoff_surf_stdev .lt. 0d0) &
   call warning(whoami,"Cutoff for density smoothing"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef
!#undef __COMPUTE
! GAUSSIAN standard deviations for smoothing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (.not.existtag_nocase("stdev_charge")) then
   call warning(whoami, "Gaussian stdev for charge smoothing"//' unspecified. Will compute from grid.',0)
   stdev_charge=-1
  else
   keyword=getval("stdev_charge");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
                        "Gaussian stdev for charge smoothing"//' to '//keyword(1:l)//'.');
   stdev_charge = atof(keyword)
   if ( stdev_charge .lt. 0d0) &
   call warning(whoami,"Gaussian stdev for charge smoothing"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef __DEFAULT
!#undef
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (.not.existtag_nocase("stdev_kappa")) then
   call warning(whoami, "Gaussian stdev for ionic strength smoothing"//' unspecified. Will compute from grid.',0)
   stdev_kappa=-1
  else
   keyword=getval("stdev_kappa");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
                        "Gaussian stdev for ionic strength smoothing"//' to '//keyword(1:l)//'.');
   stdev_kappa = atof(keyword)
   if ( stdev_kappa .lt. 0d0) &
   call warning(whoami,"Gaussian stdev for ionic strength smoothing"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef __DEFAULT
!#undef
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (.not.existtag_nocase("stdev_eps")) then
   call warning(whoami, "Gaussian stdev for dielectric smoothing"//' unspecified. Will compute from grid.',0)
   stdev_eps=-1
  else
   keyword=getval("stdev_eps");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
                        "Gaussian stdev for dielectric smoothing"//' to '//keyword(1:l)//'.');
   stdev_eps = atof(keyword)
   if ( stdev_eps .lt. 0d0) &
   call warning(whoami,"Gaussian stdev for dielectric smoothing"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef __DEFAULT
!#undef
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (.not.existtag_nocase("stdev_surf")) then
   call warning(whoami, "Gaussian stdev for density smoothing"//' unspecified. Will compute from grid.',0)
   stdev_surf=-1
  else
   keyword=getval("stdev_surf");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
                        "Gaussian stdev for density smoothing"//' to '//keyword(1:l)//'.');
   stdev_surf = atof(keyword)
   if ( stdev_surf .lt. 0d0) &
   call warning(whoami,"Gaussian stdev for density smoothing"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef __DEFAULT
!#undef
!
!%%%%%%%%%%%%%%%%%%%%% compute normalized parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 oo_stdev_charge=1d0/stdev_charge; oo_stdev_eps=1d0/stdev_eps; oo_stdev_kappa=1d0/stdev_kappa; oo_stdev_surf=1d0/stdev_surf
! oo_stdev_charge2=oo_stdev_charge**2;oo_stdev_eps2=oo_stdev_eps**2;oo_stdev_surf2=oo_stdev_surf**2;oo_stdev_kappa2=oo_stdev_kappa**2
 cutoff_eps =cutoff_eps_stdev*stdev_eps
 cutoff_kappa =cutoff_kappa_stdev*stdev_kappa
 cutoff_surf =cutoff_surf_stdev*stdev_surf
 cutoff_charge=cutoff_charge_stdev*stdev_charge
!%%%%%%%%%%%%%%%%%%%%%%%%%%% done with input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
!
  end subroutine molecule_read_parameters
!
!**********************************************************************************
  subroutine molecule_grid_objects(xcen, ycen, zcen, dxcor, dycor, dzcor, eps, kappa, rhs, nx, ny, nz)
!DEC$ ATTRIBUTES FORCEINLINE :: erfo7
!DEC$ ATTRIBUTES FORCEINLINE :: erfo5
! to avoid circular dependencies, need to compile 'molecule' before 'grid', so pass the gridded data explicitly
! use gridsize, only: nx, ny, nz
! use grid, only : xcen, ycen, zcen, dxcor, dycor, dzcor
! use state, only: eps, kappa, rhs
  implicit none
! declare passed parameters (gridded data) -- see above comments
  integer :: nx, ny, nz
  real*8, dimension(nx,ny,nz) :: eps, kappa, rhs
  real*8, dimension(nx) :: xcen, dxcor
  real*8, dimension(ny) :: ycen, dycor
  real*8, dimension(nz) :: zcen, dzcor
! local vars
  character(len=22), parameter :: whoami='MOLECULE_GRID_OBJECTS'
  integer :: l
  integer :: i, j, k, ip1, jp1, kp1, im1, jm1, km1, idir, jdir, kdir
  real*8 :: xpt, ypt, zpt,&
           dz2, dzdy2, dzdydx2,&
           d, d2, e, e2, val, qn
!
  real*8 :: reps , repsn, r2epsn , rmax2eps
  real*8 :: rkappa , rkappan, r2kappan , rmax2kappa
  real*8 :: rcharge, rchargen, r2chargen, rmax2charge
  real*8 :: rsurf , rsurfn, r2surfn , rmax2surf
!
  logical :: qkappai, qkappaj, qkappak,&
          qchargei, qchargej, qchargek,&
          qepsi, qepsj, qepsk,&
          qsurfi, qsurfj, qsurfk,&
          qanyi, qanyj, qanyk
!
  real*8 :: erfsun ! external erf function; code by Sun Microsystems
  real*8 :: erfo7 ! approximation to within O(-7)
  real*8 :: erfo5 ! approximation to within O(-5) [contains an exponential term]
!
  if (.not.molecule_initialized) call molecule_initialize()
!
! if stdevs are negative, compute from grid :
  d=-1
!
  if (oo_stdev_eps.lt.0) then;
   if (d.lt.0) d=1.5d0*min(minval(dxcor),minval(dycor),minval(dzcor))
   stdev_eps=d; cutoff_eps=cutoff_eps_stdev*stdev_eps; oo_stdev_eps=1d0/stdev_eps
  endif
  if (oo_stdev_kappa.lt.0) then;
   if (d.lt.0) d=1.5d0*min(minval(dxcor),minval(dycor),minval(dzcor))
   stdev_kappa=d; cutoff_kappa=cutoff_kappa_stdev*stdev_kappa; oo_stdev_kappa=1d0/stdev_kappa
  endif
  if (oo_stdev_charge.lt.0) then;
   if (d.lt.0) d=1.5d0*min(minval(dxcor),minval(dycor),minval(dzcor))
   stdev_charge=d; cutoff_charge=cutoff_charge_stdev*stdev_charge; oo_stdev_charge=1d0/stdev_charge
  endif
  if (oo_stdev_surf.lt.0) then;
   if (d.lt.0) d=1.5d0*min(minval(dxcor),minval(dycor),minval(dzcor))
   stdev_surf=d; cutoff_surf=cutoff_surf_stdev*stdev_surf; oo_stdev_surf=1d0/stdev_surf
  endif
! write(0,*) 'STDEV:', oo_stdev_eps, oo_stdev_surf ! aa
!
  if (.not.associated(rho)) allocate(rho(nx,ny,nz))
! initialize all but rhs (so that we can superpose additional sources outside of this module)
  do k=1,nz; do j=1, ny; do i=1,nx ;
                                       eps(i,j,k)=0d0; kappa(i,j,k)=0d0; rho(i,j,k)=0d0;
  enddo; enddo; enddo;
! gridding
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 1) for each atom, find the closest gridpoint
  allocate(closest_pt(3,natom));
! begin search in the middle (in the future, can start alternatively from existing closest_pt array)
  i=nx/2; j=ny/2; k=nz/2;
!
! process grid centers (may need to do the same for corner grid to compute gradients (which will defined naturally on the corners)
!
  do l=1, natom
!
   xpt=r(1,l)
   ypt=r(2,l)
   zpt=r(3,l)
! write(666,*) xpt,ypt,zpt
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% x-
   if ( xpt.gt.xcen(i) ) then
    do
     ip1=i+1
     if (ip1.eq.nx) then
      call warning(whoami, 'MOLECULE INTERSECTS RIGHT X-GRID BOUNDARY.',-1)
      closest_pt(1,l)=i
      exit
     elseif (xpt.le.xcen(ip1)) then
      closest_pt(1,l)=i
      exit
     endif
     i=ip1
    enddo
!
   else ! xpt.le.xcen(i)
    do
     im1=i-1
     if (im1.eq.1) then
      call warning(whoami, 'MOLECULE INTERSECTS LEFT X-GRID BOUNDARY.',-1)
      closest_pt(1,l)=im1
      exit
     elseif (xpt.ge.xcen(im1)) then
      closest_pt(1,l)=im1
      exit
     endif
     i=im1
    enddo
   endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% y-
   if ( ypt.gt.ycen(j) ) then
    do
     jp1=j+1
     if (jp1.eq.ny) then
      call warning(whoami, 'MOLECULE INTERSECTS RIGHT Y-GRID BOUNDARY.',0)
      closest_pt(2,l)=j
      exit
     elseif (ypt.le.ycen(jp1)) then
      closest_pt(2,l)=j
      exit
     endif
     j=jp1
    enddo
!
   else ! ypt.le.ycen(j)
    do
     jm1=j-1
     if (jm1.eq.1) then
      call warning(whoami, 'MOLECULE INTERSECTS LEFT Y-GRID BOUNDARY.',0)
      closest_pt(2,l)=jm1
      exit
     elseif (ypt.ge.ycen(jm1)) then
      closest_pt(2,l)=jm1
      exit
     endif
     j=jm1
    enddo
   endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% z-
   if ( zpt.gt.zcen(k) ) then
    do
     kp1=k+1
     if (kp1.eq.nz) then
      call warning(whoami, 'MOLECULE INTERSECTS RIGHT Z-GRID BOUNDARY.',0)
      closest_pt(3,l)=k
      exit
     elseif (zpt.le.zcen(kp1)) then
      closest_pt(3,l)=k
      exit
     endif
     k=kp1
    enddo
!
   else ! zpt.le.zcen(k)
    do
     km1=k-1
     if (km1.eq.1) then
      call warning(whoami, 'MOLECULE INTERSECTS LEFT Z-GRID BOUNDARY.',0)
      closest_pt(3,l)=km1
      exit
     elseif (zpt.ge.zcen(km1)) then
      closest_pt(3,l)=km1
      exit
     endif
     k=km1
    enddo
   endif
!
  enddo ! over all atoms
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! loop over all points in the vicinity of the closest point (defined by distance)
  do l=1, natom ! this loop should be parallelized
   xpt=r(1,l)
   ypt=r(2,l)
   zpt=r(3,l)
!
! write(777,*) radius(l), q(l) ! aa
! write(888,*) r(:,l), radius(l), closest_pt(:,l) ! aa; coors correct
! precomputation of parameters for the smoothed densities
   reps= (radius(l)+eps_padding ); repsn =reps * oo_stdev_eps ! normalize by stdev
   rkappa= (radius(l)+kappa_padding ); rkappan =rkappa * oo_stdev_kappa
   rcharge=(radius(l)+charge_padding); rchargen=rcharge * oo_stdev_charge
   rsurf =(radius(l)+surf_padding ); rsurfn =rsurf * oo_stdev_surf
!
! r2epsn =repsn **2
! r2kappan =rkappan **2
! r2chargen=rchargen**2
! r2surfn =rsurfn **2
!
   rmax2eps =(cutoff_eps +reps )**2
   rmax2kappa =(cutoff_kappa +rkappa )**2
   rmax2charge=(cutoff_charge+rcharge)**2
   rmax2surf =(cutoff_surf +rsurf )**2
! compute nominal charge density ( nominal because it will not be accurate on the grid and will need correction )
   qn=q(l) * 0.75d0 / ( pi * rchargen **3 )
!
   do kdir=-1,1,2 ! directions
    k=closest_pt(3,l)+(kdir+1)/2; ! when kdir = 1, shift to the the right relative to closest point, but not when kdir=-1
    qepsk=.true.; qkappak=.true.; qchargek=.true.; qsurfk=.true.
    do
! z-grid OOB test:
     if ((k.ge.nz).or.(k.le.1)) then
      call warning(whoami, 'SMOOTHED MOLECULE INTERSECTS Z-GRID BOUNDARY.',0)
      exit
     endif
!
     dz2=(zcen(k)-zpt)**2
! k range test:
     if (qepsk) qepsk=(dz2.lt.rmax2eps)
     if (qkappak) qkappak=(dz2.lt.rmax2kappa)
     if (qchargek) qchargek=(dz2.lt.rmax2charge)
     if (qsurfk) qsurfk=(dz2.lt.rmax2surf )
     qanyk = qepsk .or. qkappak .or. qchargek .or. qsurfk
     if (.not.( qanyk )) exit
!
     do jdir=-1,1,2
      j=closest_pt(2,l)+(jdir+1)/2;
      qepsj=qepsk; qkappaj=qkappak; qchargej=qchargek; qsurfj=qsurfk
      do
! y-grid OOB test:
       if ((j.ge.ny).or.(j.le.1)) then
        call warning(whoami, 'SMOOTHED MOLECULE INTERSECTS Y-GRID BOUNDARY.',0)
        exit
       endif
!
       dzdy2=dz2+(ycen(j)-ypt)**2
! j range test:
       if (qepsj) qepsj=(dzdy2.lt.rmax2eps)
       if (qkappaj) qkappaj=(dzdy2.lt.rmax2kappa)
       if (qchargej) qchargej=(dzdy2.lt.rmax2charge)
       if (qsurfj) qsurfj=(dzdy2.lt.rmax2surf )
       qanyj = qepsj .or. qkappaj .or. qchargej .or. qsurfj
       if (.not.( qanyj )) exit
!
       do idir=-1,1,2
        i=closest_pt(1,l)+(idir+1)/2;
        qepsi=qepsj; qkappai=qkappaj; qchargei=qchargej; qsurfi=qsurfj
        do
! x-grid OOB test:
         if ((i.ge.nx).or.(i.le.1)) then
          call warning(whoami, 'SMOOTHED MOLECULE INTERSECTS X-GRID BOUNDARY.',0)
          exit
         endif
!
         dzdydx2=dzdy2+(xcen(i)-xpt)**2
! i range test:
         if (qepsi) qepsi=(dzdydx2.lt.rmax2eps)
         if (qkappai) qkappai=(dzdydx2.lt.rmax2kappa)
         if (qchargei) qchargei=(dzdydx2.lt.rmax2charge)
         if (qsurfi) qsurfi=(dzdydx2.lt.rmax2surf )
         qanyi = qepsi .or. qkappai .or. qchargei .or. qsurfi
         if (.not.( qanyi )) exit
         d=sqrt(dzdydx2)
!
!#undef 1
! grid epsilon %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if (qepsi &
                     .and. ( eps(i,j,k) .lt. 1.2 ) &
         ) then ! do not compute if cutoff exceeded in previous iteration
! compute atomic contribution to mesh point
           d2=d*oo_stdev_eps ; e = osq2 * (d2 + repsn) ; e2= osq2 * (d2 - repsn)
!
           val=half * ( erfo7 (e) - erfo7 (e2) ) + osq2pi*( exp(-e**2) - exp(-e2**2 ) ) /d2
!
           eps(i,j,k) = eps(i,j,k) + val ! scale to correct levels later
!
! this is the convolution of a normalized gaussian with stdev sigma with a function that changes across the boundary of
! a sphere of the given atomic radius (function values depend on the property: charge, surf, kappa, eps)
         endif ! epsilon
! grid kappa %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if (qkappai &
                     .and. ( kappa(i,j,k) .lt. 1.2 ) &
         ) then ! do not compute if cutoff exceeded in previous iteration
! compute atomic contribution to mesh point
           d2=d*oo_stdev_kappa ; e = osq2 * (d2 + rkappan) ; e2= osq2 * (d2 - rkappan)
!
           val=half * ( erfo7 (e) - erfo7 (e2) ) + osq2pi*( exp(-e**2) - exp(-e2**2 ) ) /d2
!
           kappa(i,j,k)=kappa(i,j,k) + val
         endif ! kappa
! grid charges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if (qchargei) then ! do not compute if cutoff exceeded in previous iteration
! compute atomic contribution to mesh point
           d2=d*oo_stdev_charge ; e = osq2 * (d2 + rchargen) ; e2= osq2 * (d2 - rchargen)
!
           val=half * ( erfo7 (e) - erfo7 (e2) ) + osq2pi*( exp(-e**2) - exp(-e2**2 ) ) /d2
!
           rhs(i,j,k)=rhs(i,j,k) + val * qn ! qn is nominal charge density: q / ( 4/3pi r^3 )
         endif ! charge
! grid surface %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if (qsurfi &
! this is probably not a good idea if the surface derivatives are needed
!#ifdef 1
! .and. ( rho(i,j,k) .lt. 1.2 ) &
!#endif
         ) then ! do not compute if cutoff exceeded in previous iteration
! compute atomic contribution to mesh point
           d2=d*oo_stdev_surf ; e = osq2 * (d2 + rsurfn) ; e2= osq2 * (d2 - rsurfn)
!
           val=half * ( erfo7 (e) - erfo7 (e2) ) + osq2pi*( exp(-e**2) - exp(-e2**2 ) ) /d2
!
           rho(i,j,k)=rho(i,j,k) + val ! density: 1 inside sphere, 0 outside sphere
!
         endif ! surface
!
         i=i+idir
        enddo ! over i
       enddo ! idir changes from -1 to 1
!
       j=j+jdir
      enddo
     enddo ! jdir changes from -1 to 1
!
     k=k+kdir
    enddo
   enddo ! kdir changes from -1 to 1
!if (l.eq.100) exit ! aa -- grid a few atoms
  enddo ! over all atoms
!
! 3) Need to scale gridded coefficients so that they add up correctly
! The most straightforward and, perhaps, correct way to do this is
! to keep track of the grid points that support a given charge, and
! scale all of them uniformly so that the charge density integrated over this support
! integrates to the charge.
! Another approach, which is adopted below, is to scale all of the
! grid values uniformly so that the charge density integrated over the entire domain
! adds up to the total charge
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! scale charges uniformly
  qn=0d0 ! total grid charge
! scaling factors for kappa and eps
  d =eps_solute-eps_solvent
  d2=kappa_solute-kappa_solvent
!
  do k=2, nz-1
   dz2=dzcor(k)
   do j=2, ny-1
    dzdy2=dz2*dycor(j)
    do i=2, nx-1
!
     qn = qn + rhs(i,j,k) * dxcor(i) * dzdy2
!
! scale rho, eps and kappa
!%%%%%%%%%%%%%%%%%%%%%%%%%
! surface
     if (rho(i,j,k).gt.1d0) rho(i,j,k)=1d0
! dielectric
     if (eps(i,j,k).gt.1d0) eps(i,j,k)=1d0 ; eps (i,j,k)=eps (i,j,k) * d + eps_solvent
! ionic strength
     if (kappa(i,j,k).gt.1d0) kappa(i,j,k)=1d0 ; kappa(i,j,k)=kappa(i,j,k) * d2 + kappa_solvent
!
!%%%%%%%%%%%%%%%%%%%%%%%%%
    enddo
    eps(1,j,k) = eps(1,j,k)*d + eps_solvent; eps(nx,j,k) = eps(nx,j,k)*d + eps_solvent
    kappa(1,j,k) = kappa(1,j,k)*d + kappa_solvent; kappa(nx,j,k) = kappa(nx,j,k)*d + kappa_solvent
   enddo
!
   do i=1,nx
    eps(i,1,k) = eps(i,1,k)*d + eps_solvent; eps(i,ny,k) = eps(i,ny,k)*d + eps_solvent
    kappa(i,1,k) = kappa(i,1,k)*d + kappa_solvent; kappa(i,ny,k) = kappa(i,ny,k)*d + kappa_solvent
   enddo
!
  enddo
!
  do j=1,ny ; do i=1,nx
   eps(i,j,1) = eps(i,j,1)*d + eps_solvent; eps(i,j,nz) = eps(i,j,nz)*d + eps_solvent
   kappa(i,j,1) = kappa(i,j,1)*d + kappa_solvent; kappa(i,j,nz) = kappa(i,j,nz)*d + kappa_solvent
  enddo ; enddo
!
  rhs=rhs*sum(q)/qn;
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! done !
  end subroutine molecule_grid_objects
!**********************************************************************************
  subroutine molecule_initialize()
  use PDB
  implicit none
!
  character(len=19) , parameter :: whoami='MOLECULE_INITIALIZE'
!
  character(len=20) :: struct_type, charge_type, radius_type, coords_type
  character(len=10) :: column ! specifies from which column (occupancy or B) charge/radius data are read
  character(len=10) :: keyword
  character(len=100) :: parmfilename, structfilename, coorfilename, chargefilename, radiusfilename
!
  logical :: need_coordinates=.true., need_charges=.true., need_radii=.true.
!
  integer :: i
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (molecule_initialized) call molecule_done()
!
  struct_type=getval_nocase('STRUCTURE_FORMAT'); call toupper(struct_type);
  structfilename=getval_nocase('STRUCTURE'); ! structure file
!
  select case(struct_type);
   case('PSF');
    call system_read_structure(structfilename); ! in dynamol
    need_coordinates=.true.
    need_charges=.false. ! charges are specified in PSF, but can be specified optionally
    need_radii=.true. ! charges are specified in PSF, but can be specified optionally
! see also if parameter file(s) have been specified
    if (existtag_nocase('PARAMETERS')) then
!%%%%%%%%%%%%%%%%%%% read parameter file(s)%%%%%%%%%%%%%%%%%%%%%%%%
     parmfilename=getval_nocase('PARAMETERS') ! parameter file(s)
     call system_read_parameters(parmfilename)
     i=2
     do
      write(keyword,'(I10)') i
      call adjustleft(keyword)
      if (existtag_nocase('PARAMETERS'//keyword(1:len_trim(keyword)))) then
        parmfilename=getval_nocase('PARAMETERS'//keyword(1:len_trim(keyword)))
        call system_read_parameters(parmfilename)
        i=i+1
      else
        exit
      endif
     enddo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! DISABLE CODE BELOW: COMPUTING RADII REQUIRES THAT COORDINATES BE INITIALIZED WHICH IS FALSE AT THIS STAGE
! if (system_parameters_initialized) then
! call system_get_vw_radius()
! need_radii=.false.
! endif
    endif
! other structure files
   case('PDB');
    call PDB_init(structfilename,'PDB'); ! last argument specifies whether PDB or PQR format will be used
    need_coordinates=.false.
    need_radii=.true.
    need_charges=.true.
   case('PQR');
    call PDB_init(structfilename,'PQR');
    need_coordinates=.false.
    need_radii=.false.
    need_charges=.false.
   case('CHARMM');
    call PDB_init(structfilename,'CHARMM');
    need_coordinates=.false.
    need_radii=.true.
    need_charges=.true.
   case default
    call warning(whoami, 'Unknown structure format. Abort',-1)
    return
  end select
! coordinates
  if (need_coordinates.or.existtag_nocase('COORDINATES')) then
   coorfilename=getval_nocase('COORDINATES');
   call system_read_coordinates(coorfilename); ! in dynamol; this will work because PDB_init hacks the system module not to die in the absence of
                                               ! proper topology and parameters
  elseif (need_coordinates.and..not.existtag_nocase('COORDINATES')) then
   call error(whoami, 'Atomic coordinates unspecified. Abort',-1);
  endif
!
! charges
  if (need_charges.or.existtag_nocase('charges')) then
   charge_type=getval_nocase('chargestype'); call toupper(charge_type);
   chargefilename=getval_nocase('charges')
!
   select case(charge_type);
    case('PDB');
     column=getval_nocase('charges_column');
     call PDB_read_charges(chargefilename,'PDB',column)
    case('PQR');
     call PDB_read_charges(chargefilename,'PQR')
    case('ATOMID');
     call PDB_read_charges(chargefilename,'ATOMID')
    case('CHARMM');
     call PDB_read_charges(chargefilename,'CHARMM') ! charges are assumed to be in main weighting array
    case default
     call error(whoami, 'Charges in unknown format. Abort',-1)
    return
   end select
  elseif (need_charges.and..not.existtag('charges')) then
   call error(whoami, 'Atomic charges unspecified. ABORT',-1);
  endif
!
! radii
  if (need_radii.or.existtag('radii')) then
   radius_type=getval('radiitype'); call toupper(radius_type);
!
   select case(radius_type);
    case('PDB');
     column=getval('radii_column');
     radiusfilename=getval('radii')
     call PDB_read_radii(radiusfilename,'PDB',column)
    case('PQR');
     radiusfilename=getval('radii')
     call PDB_read_radii(radiusfilename,'PQR')
    case('ATOMID');
     radiusfilename=getval('radii')
     call PDB_read_radii(radiusfilename,'ATOMID')
    case('CHARMM');
     radiusfilename=getval('radii')
     call PDB_read_radii(radiusfilename,'CHARMM') ! radii are assumed to be in main weighting array
    case('PARAM','PARAMETER','PARM'); ! obtain from parameter file; this works only if proper structure file is present:
     if (struct_type.eq.'PSF') then
      if(.not.system_parameters_initialized) then ! if params read above do nothing
!%%%%%%%%%%%%% read parameter file(s) - note code duplication above %
       parmfilename=getval_nocase('parameters') ! parameter file(s)
       call system_read_parameters(parmfilename)
       i=2
       do
        write(keyword,'(I10)') i
        call adjustleft(keyword)
        if (existtag_nocase('parameters'//keyword(1:len_trim(keyword)))) then
          parmfilename=getval_nocase('parameters'//keyword(1:len_trim(keyword)))
          call system_read_parameters(parmfilename)
          i=i+1
        else
          exit
        endif
       enddo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      endif
      call system_get_vw_radius()
     else ! not PSF
      call error(whoami, 'CANNOT OBTAIN RADII FROM PARAMETER FILE(S) BECAUSE ATOM TYPES ARE MISSING FROM STRUCTURE. ABORT.',-1)
     endif
    case default
     call error(whoami, 'RADII IN UNKNOWN FORMAT. ABORT.',-1)
    return
   end select
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif (need_radii.and..not.existtag('RADII')) then
   call error(whoami, 'ATOMIC RADII UNSPECIFIED. ABORT',-1);
  endif
!
  nullify(rho)
  molecule_initialized=.true.
!
  end subroutine molecule_initialize
!***********************************************************************************
  subroutine molecule_center(rcen, qmass)
  use sysmanip
  use sysinfo
  implicit none
  real*8 :: dimens(ndim,3)
  real*8 :: rcen(3)
  logical :: qmass
!
  if (.not. molecule_initialized) call molecule_initialize()
!
  dimens=sysinfo_dimens(qmass,(/-1/)) ! last argument is an array of length 1
  call sysmanip_translate(rcen-dimens(1:ndim,1),(/-1/))
!
  end subroutine molecule_center
!***********************************************************************************
  subroutine molecule_align(qmass)
  use sysmanip
  use sysinfo
  implicit none
  logical :: qmass
!
  if (.not. molecule_initialized) call molecule_initialize()
!
  call sysmanip_align_pc(qmass,(/-1/),(/-1/))! align the principal components of the molecule with the Cartesian vectors
!
  end subroutine molecule_align
!***********************************************************************************
  function molecule_dimens()
  use sysinfo
  implicit none
  real*8 :: molecule_dimens(ndim,3)
  real*8 :: pad
!
  if (.not. molecule_initialized) call molecule_initialize()
!
  molecule_dimens=sysinfo_dimens(.false.,(/-1/)) ! last argument is an array of length 1
! increase maximum coordinate values computed above by the padding
  pad=max(eps_padding+cutoff_eps, kappa_padding+cutoff_kappa, charge_padding+cutoff_charge, surf_padding+cutoff_surf)
  molecule_dimens(:,2)=molecule_dimens(:,2)-pad;
  molecule_dimens(:,3)=molecule_dimens(:,3)+pad;
!
  end function molecule_dimens
!
!***********************************************************************************
  integer function molecule_ndim() ; implicit none
  molecule_ndim=ndim
  end function molecule_ndim
!***********************************************************************************
  subroutine molecule_done()
  implicit none
  molecule_initialized=.false.
  if (allocated(closest_pt)) deallocate(closest_pt)
  if (associated(rho)) deallocate(rho)
  call system_done()
  end subroutine molecule_done
!***********************************************************************************
end module molecule
