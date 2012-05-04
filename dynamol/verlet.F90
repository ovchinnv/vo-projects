! **********************************************************************!
! This source file was was generated automatically from a master source !
! code tree, which may or may not be distributed with this code, !
! because it is up to the distributor, and not up to me. !
! If you edit this file (rather than the master source file) !
! your changes will be lost if another pull from the master tree occurs.!
! In case you are wondering why, this approach makes it possible for !
! me to have the same master source code interfaced with different !
! applications (some of which are written in a way that is quite far !
! from being object-oriented) at the source level. !
! **********************************************************************!
module verlet
 private
 logical, save :: langevin ! indicates that Langevin dynamics active
 integer, save :: random_par=1 ! channel from which to draw Random numbers; must not be read-only (see clcg.f)
!
 real*8, save :: dtau ! timestep (fs)
!
! variables for Verlet & Langevin dynamics
 real*8, pointer, dimension(:), save :: langevin_friction, langevin_temperature
! LD integration variables
 real*8, allocatable, dimension(:,:), save :: noise, halfdtom, rho, a, b
 logical :: verlet_initialized=.false.
 public verlet_init
 public verlet_integrate
 public verlet_done
 contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine verlet_init()
  use parser
  use output
  use rng
  use constants
!
  use system
!
  implicit none
!
 character(len=11) :: whoami='VERLET_INIT'
 integer :: i
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (verlet_initialized) call verlet_done()
  nullify(langevin_temperature, langevin_friction)
!
  if (.not.(system_coordinates_initialized.and.&
             system_velocities_initialized.and.&
             system_ok)) then
    call error(whoami,' SYSTEM NOT INITIALIZED',-1)
    return
  endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  dtau=atof(getval('dt')) ! timestep
  if (dtau.lt.0d0) then
    call error(whoami, ' TIMESTEP IS NEGATIVE',-1)
    return
  else
    dtau=dtau*akma_per_fs
  endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  allocate(halfdtom(ndim,natom), a(ndim,natom), b(ndim,natom))
!
  langevin=atol(getval_nocase('langevin')) ! is langevin dynamics active?
  if (langevin) then
   langevin_temperature=>system_init_mol_scalar_nocase('langevin_temperature')
   if (.not.associated(langevin_temperature)) then
    call warning(whoami, 'ERROR INITIALIZING "LANGEVIN_TEMPERATURE"',-1)
   elseif (any(langevin_temperature.lt.zero)) then
    call warning(whoami, 'LANGEVIN TEMPERATURE ARRAY CONTAINS NEGATIVE VALUE(S)',-1)
   else
    langevin_friction=>system_init_mol_scalar_nocase('langevin_friction')
    if (.not.associated(langevin_friction)) then
     call warning(whoami, 'ERROR INITIALIZING "LANGEVIN_FRICTION"',-1)
    elseif (any(langevin_friction.lt.zero)) then
     call warning(whoami, 'LANGEVIN FRICTION ARRAY CONTAINS NEGATIVE VALUE(S)',-1)
    else
     langevin_friction=langevin_friction*0.001d0*fs_per_akma ! input assumed to be in ps^(-1)
     allocate(noise(ndim,natom), rho(ndim,natom))
! precompute random numbers for first verlet step
     call randomg_vector(noise,natom*ndim,random_par) ! assume RGN has been initialized
    endif ! friction
   endif ! temperature
  else
   allocate(langevin_friction(natom)); langevin_friction=zero
  endif ! langevin
! quit if the warnings above are fatal
! termination will be invoked at a different level
!
  if (.not.fatal_warning()) then
!
   do i=1,ndim
    a(i,:)=-half*langevin_friction*dtau+one;
    b(i,:)=one/(-a(i,:)+two);
    halfdtom(i,:)=half*dtau/m
    if (langevin) rho(i,:)=sqrt(langevin_friction*langevin_temperature*kboltzmann/halfdtom(i,:));
   enddo
!
   call system_compute()
!
   verlet_initialized=.true.
!
  else
   deallocate(a,b,halfdtom)
  endif
! have no need for the temperature or friction arrays
!
  if (associated(langevin_friction)) deallocate(langevin_friction)
  if (associated(langevin_temperature)) deallocate(langevin_temperature)
!
!
  end subroutine verlet_init
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine verlet_done()
   implicit none
   if (associated(langevin_temperature)) deallocate(langevin_temperature)
   if (associated(langevin_friction)) deallocate(langevin_friction)
!
   if (allocated(halfdtom)) deallocate(halfdtom)
   if (allocated(rho)) deallocate(rho)
   if (allocated(a)) deallocate(a)
   if (allocated(b)) deallocate(b)
   if (allocated(noise)) deallocate(noise)
   verlet_initialized=.false.
  end subroutine verlet_done
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine verlet_integrate(nsteps)
! velocity verlet
! not sure what the reference is; this seems to be similar to BBK
! we are solving:
!
! "                   .
! m x = - grad(V) - g m x + SQRT[ 2 g k T m d] eta
! eta(t) is Gaussian white noise with zero mean and unit variance
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   use rng
   use system
!
   implicit none
   integer :: j, nsteps
!
   if (.not.verlet_initialized) call verlet_init()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   do j=1,nsteps
!%%%%%%%%%%%%%%%%%%% half a kick
    if (langevin) then
      vr=a*vr+halfdtom*(fr+rho*noise)
    else
      vr=a*vr+halfdtom*fr
    endif ! langevin
!%%%%%%%%%%%%%%%%%%%% drift
    r=r+dtau*vr
!% compute new forces and noise
    call system_compute()
!%%%%%%%%%%%%%%%%%%%% half a kick
    if (langevin) then
      call randomg_vector(noise,natom*ndim,random_par)
      vr=b*(vr+halfdtom*(fr+rho*noise))
    else ! langevin
      vr=b*(vr+halfdtom*fr)
    endif ! langevin
!
   enddo ! n
!
  end subroutine verlet_integrate

end module verlet
