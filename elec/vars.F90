!
!
 module vars
 public
!
! __IPAR(nx,16)
! __IPAR(ny,16)
! __IPAR(nz,16)
!
 __IPAR(nx,32)
 __IPAR(ny,32)
 __IPAR(nz,32)
!
! __IPAR(nx,48)
! __IPAR(ny,48)
! __IPAR(nz,48)
!
! __IPAR(nx,64)
! __IPAR(ny,64)
! __IPAR(nz,64)

! domain size
 __FPAR(x0,-1d0)
 __FPAR(y0,-1d0)
 __FPAR(z0,-1d0)
 __FPAR(x1,1d0)
 __FPAR(y1,1d0)
 __FPAR(z1,1d0)
 __FPAR(Lx,x1-x0)
 __FPAR(Ly,y1-y0)
 __FPAR(Lz,z1-z0)
 __FPAR(oLx,1d0/Lx)
 __FPAR(oLy,1d0/Ly)
 __FPAR(oLz,1d0/Lz)
 __FPAR(oV, oLx*oLy*oLz)

! grid
 __FARR(xx,:)
 __FARR(yy,:)
 __FARR(zz,:)

! positions and charges
 int :: npt=0
 __FARR(x,:);
 __FARR(y,:);
 __FARR(z,:);
 __FARR(q,:);
! electrostatic potential variables
 __FARR(rho,:,:,:);   ! smoothed density
 __FARR(phi,:,:,:);   ! L/R potential by ewald summation
 __FARR(lap,:,:,:);   ! laplacian for Finite Difference solution
 __FARR(phisr_q,:);   ! short range potential and gradients below
 __FARR(phisr_dx_q,:);
 __FARR(phisr_dy_q,:);
 __FARR(phisr_dz_q,:);
 
 float, target :: ene(3)=0d0
 float, pointer :: el=>ene(1) ! long-range electrostatic potential
 float, pointer :: elsr=>ene(2) ! short-range electrostatic potential
 float, pointer :: el_self=>ene(3) ! self-interaction part of long-range energy
!
 __FARR(grad_el,:,:)
 __FARR(el_dx,:) ; ! gradient of l/r electrostatic potential; will be bound to grad_el
 __FARR(el_dy,:) ; !
 __FARR(el_dz,:) ; !
!
 __FARR(grad_elsr,:,:)
 __FARR(elsr_dx,:) ; ! gradient of s/r electrostatic potential; will be bound to grad_elsr
 __FARR(elsr_dy,:) ; !
 __FARR(elsr_dz,:) ; !
!
 __IPAR(strmax,1000);
 character(len=strmax) :: chargefile='xyzq.dat'
 int :: chunit=100

 __FPAR(pi,3.141592653589793d0)
 __FPAR(twopi,2d0*pi)
 __FPAR(fourpi,4d0*pi)
! should be softcoded below 
 __FPAR(eps,1d0); ! permittivity
 __FPAR(oeps,1d0/eps); ! permittivity
 __FPAR(spt,1.25d0); ! filter support -- should move this elsewhere and make soft
 __FPAR(oos,1d0/spt);

! for long-range ewald
  __IPAR(nwave,20) ! wavenumber range
#if _FILTER==_GAUSS
 __FPAR(ktol,1d-5)
#else
 __FPAR(ktol, 0.5d0*0.9d0**fpo) ! empirical error tolerance for when to switch to taylor expansion
#endif
! __FPAR(ktol, 0.2d0) ! empirical error tolerance for when to switch to taylor expansion
 __FPAR(kzero, 1d-14);
!
 bool :: quiet=.false. ! flag to generate minimal output

 end module vars
