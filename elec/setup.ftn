
 subroutine setup
  use vars
  int :: iostatus=0
  int i, lines

  open(file=chargefile, unit=chunit, form='formatted', status='old', iostat=iostatus)
  if (iostatus.ne.0) then
   __ERR('Could not read file ', trim(chargefile))
   __DIE
  endif
!
! read to the end of file to determine number of particles
!
  lines=0
  do while (iostatus==0)
   read(chunit,*,iostat=iostatus)
   __INC(lines)
  enddo
  __INCR(lines,-1)
  __OUT('file "', trim(chargefile), '" contains ',lines,'lines' )
! now read data
  allocate(x(lines), y(lines), z(lines), q(lines))
  rewind(chunit)
  do i=1, lines
   read(chunit,*,iostat=iostatus) x(i), y(i), z(i), q(i)
  enddo
  npt=lines
  __OUT('read ', npt, ' coordinates and charges')
! perturb charge for FD ; this is a way to valudate del_dx, etc.
! x(2)=x(2)+0.001;

  __OUT('grid dimensions are ', nx, ny, nz)
  allocate(xx(nx), yy(ny), zz(nz)) ! predetermined grid sizes
! compute grids :
  dxx=Lx/(nx-1);
  do i=1,nx
   xx(i)=x0+(i-1)*dxx
  enddo
!
  dyy=Ly/(ny-1);
  do i=1,ny
   yy(i)=y0+(i-1)*dyy
  enddo
!
  dzz=Lz/(nz-1);
  do i=1,nz
   zz(i)=z0+(i-1)*dzz
  enddo
! write grids : 
  call write_field(xx, 'xx.txt',ascii,nx)
  call write_field(yy, 'yy.txt',ascii,ny)
  call write_field(zz, 'zz.txt',ascii,nz)
!
  __ALLOC(rho(nx,ny,nz));   ! smoothed density
  __ALLOC(phisr_q(npt));    ! short range potential and gradienst below
  __ALLOC(phisr_dx_q(npt));
  __ALLOC(phisr_dy_q(npt));
  __ALLOC(phisr_dz_q(npt));
! l/r potential
  __ALLOC(grad_el(npt,3));
  el_dx=>grad_el(:,1)
  el_dy=>grad_el(:,2)
  el_dz=>grad_el(:,3)
! s/r potential
  __ALLOC(grad_elsr(npt,3));
  elsr_dx=>grad_elsr(:,1)
  elsr_dy=>grad_elsr(:,2)
  elsr_dz=>grad_elsr(:,3)
  __ALLOC(lap(nx,ny,nz)); ! laplacian for FD Poisson

! self energy part of l/r energy can be computed right away
! note that is does not contribute to any gradients
  el_self = 0.5d0 * philr_selfc * sum(q**2) * oos * oeps;

 end subroutine setup
