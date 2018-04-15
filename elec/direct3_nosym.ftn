
subroutine direct3_nosym()
! ewald summation for long-range field split using a polynomial filter
! for debugging purposes, treat all wavenumbers separately, in particular, the negative permutations
 use vars

 int in, jn, kn, i, j
 float ik, ik2, jk, jk2, kk, kk2, kabs, kabs2
 float efac, Sr, Si, s, xk, yk, zk
 __FARR(ss,:)
 __FARR(cs,:)

 __ALLOC(ss(npt))
 __ALLOC(cs(npt))

 el=0d0 ! long range potential energy
 grad_el=0d0; ! initialize gradient

 do in=-nwave,nwave
  __OUT(' Computing contribution from wavenumber ',in,' of ', nwave);
  ik=in*twopi*oLx;
  ik2=ik**2;
  do jn=-nwave,nwave;
   jk=jn*twopi*oLy;
   jk2=jk**2;
!
   do kn=-nwave,nwave
    kk=kn*twopi*oLz;
    kabs2=kk**2+jk2+ik2;
!
    if (kabs2<kzero) cycle ! skip 0 wavenumbers
!
    kabs=sqrt(kabs2);
!
    if (kabs < ktol ) then
! use Taylor expansion
     efac = ftt(kabs)/kabs2; !
    else
! use exact FT
     efac= ft(kabs)/kabs2;
    endif
! Compute potential energy and gradient at charge points
! self-interactions to be subtracted in a different step
    Sr=0d0 ; Si=0d0;
    do i=1,npt ! precompute structure factors
     xk=ik*x(i);
     yk=jk*y(i);
     zk=kk*z(i);
! compute scalar product k . rj
     s=xk+yk+zk;
     ss(i)=sin(s);
     cs(i)=cos(s);
     __INCR(Sr, q(i) * cs(i))
     __INCR(Si, q(i) * ss(i))
    enddo ! npt
! update potential energy
    __INCR(el, efac * ( Sr**2 + Si**2 ))
! update derivatives
    do i=1,npt
! x-derivative (do not forget factor of x2 for complex conjugate, i.e. taking 2*Re( ) \)
! precompute scalar product : 
     s = cs(i)*Si - ss(i)*Sr
     __INCR( el_dx(i), efac * q(i) * ik * s );
     __INCR( el_dy(i), efac * q(i) * jk * s );
     __INCR( el_dz(i), efac * q(i) * kk * s );
    enddo ! npt
   enddo !kn
  enddo !jn
 enddo !in
! normalize properly electrostatic potential and gradients ; no x2 here because el has 1/2 prefactor ;
 el = 0.5d0 * el * oeps * oV
 el_dx = el_dx * oeps * oV ! needs 2 because : (i) account for -kx wavenumbers (2) |S(k)|^2 derivative has a 2 prefactor ; 1/2 in front of el. pot.
 el_dy = el_dy * oeps * oV
 el_dz = el_dz * oeps * oV
!
 __FREE(ss)
 __FREE(cs)
end subroutine direct3_nosym
