
subroutine direct3()
! ewald summation for long-range field split using a polynomial filter
 use vars

 int in, jn, kn, ci, cj, ck, i, j
 float ik, ik2, jk, jk2, kk, kk2, kabs, kabs2
 float efac, Sr(4), Si(4), s(4), xk, yk, zk
 __FARR(ss,:,:)
 __FARR(cs,:,:)

 __ALLOC(ss(npt,4))
 __ALLOC(cs(npt,4))

 el=0d0 ! long range potential energy
 grad_el=0d0; ! initialize gradient

 do in=0,nwave
  if (.not. quiet) __OUT(' Computing contribution from wavenumber ',in,' of ', nwave);
  ik=in*twopi*oLx;
  ik2=ik**2;
  ci=2-min(1,in) ; ! symmetry coefficients
  do jn=0,nwave;
   jk=jn*twopi*oLy;
   jk2=jk**2;
   cj=2-min(1,jn) ;
!
   do kn=0,nwave
    kk=kn*twopi*oLz;
    kabs2=kk**2+jk2+ik2;
!
    if (kabs2<kzero) cycle ! skip 0 wavenumbers
!
    kabs=sqrt(kabs2);
    ck=2-min(1,kn) ;
!
    if (kabs*spt < ktol ) then
! use Taylor expansion
     efac = ftt(kabs*spt)/kabs2/(ci*cj*ck); !
    else
! use exact FT
     efac= ft(kabs*spt)/kabs2/(ci*cj*ck);
    endif
! Compute potential energy and gradient at charge points
! self-interactions to be subtracted in a different step
    Sr=0d0 ; Si=0d0;
    do i=1,npt ! precompute structure factors
     xk=ik*x(i);
     yk=jk*y(i);
     zk=kk*z(i);
! compute scalar product k . rj
!    there are 8 wavenumbers treated here :
!                     kx ky kz
     s(1)=xk+yk+zk; ! 1  1  1
     s(2)=xk+yk-zk; ! 1  1 -1
     s(3)=xk-yk+zk; ! 1 -1  1
     s(4)=xk-yk-zk; ! 1 -1 -1
! implicitly including wavenumbers with -kx
!     s(5)=-s(4) ;  ! -1  1  1
!     s(6)=-s(3) ;  ! -1  1 -1
!     s(7)=-s(2) ;  ! -1 -1  1
!     s(8)=-s(1) ;  ! -1 -1 -1
!
! this accounts for 8 wavenumbers with the same magnitude
     do j=1,4
      ss(i,j)=sin(s(j));
      cs(i,j)=cos(s(j));
      __INCR(Sr(j), q(i) * cs(i,j))
      __INCR(Si(j), q(i) * ss(i,j))
     enddo
! implicitly including -kx wavenumbers :
!     do j=5,8
!      ss(i,j)=sin(s(j)); ! = -sin(s(8-j+1)) = -ss(i,8-j+1)
!      cs(i,j)=cos(s(j)); ! =  cos(s(8-j+1)) =  cs(i,8-j+1)
!      __INCR(Sr(j), q(i) * cs(i,j)) ! =  Sr(8-j+1)
!      __INCR(Si(j), q(i) * ss(i,j)) ! = -Si(8-j+1)
!     enddo
    enddo ! npt
! update potential energy
    __INCR(el, efac * ( sum(Sr**2) + sum(Si**2) ) ) ! omitting factor x2 for -kx wavenumbers
! update derivatives
    do i=1,npt
! x-derivative (do not forget factor of x2 for complex conjugate, i.e. taking 2*Re( ) \)
! also x2 for negative wavenumbers
! precompute scalar product : 
     do j=1,4
      s(j) = cs(i,j)*Si(j) - ss(i,j)*Sr(j)
     enddo
     __INCR( el_dx(i), efac * q(i) * ik * ( s(1) + s(2) + s(3) + s(4) ) )! &
!&           + ( -ss(i,1)*Sr(1) + cs(i,1)*Si(1) ) &
!&           + ( -ss(i,2)*Sr(2) + cs(i,2)*Si(2) ) &
!&           + ( -ss(i,3)*Sr(3) + cs(i,3)*Si(3) ) &
!&           + ( -ss(i,4)*Sr(4) + cs(i,4)*Si(4) ) &
!&         -  (  ss(i,4)*Sr(4) - cs(i,4)*Si(4))& !5
!&         -  (  ss(i,3)*Sr(3) - cs(i,3)*Si(3))& !6
!&         -  (  ss(i,2)*Sr(2) - cs(i,2)*Si(2))& !7
!&         -  (  ss(i,1)*Sr(1) - cs(i,1)*Si(1))& !8
!&       ) )
     __INCR( el_dy(i), efac * q(i) * jk * ( s(1) + s(2) - ( s(3) + s(4) ) ) )! &
!&           + ( -ss(i,1)*Sr(1) + cs(i,1)*Si(1))&
!&           + ( -ss(i,2)*Sr(2) + cs(i,2)*Si(2))&
!&           - ( -ss(i,3)*Sr(3) + cs(i,3)*Si(3))&
!&           - ( -ss(i,4)*Sr(4) + cs(i,4)*Si(4))&
!&          + (  ss(i,4)*Sr(4) - cs(i,4)*Si(4))&
!&          + (  ss(i,3)*Sr(3) - cs(i,3)*Si(3))&
!&          - (  ss(i,2)*Sr(2) - cs(i,2)*Si(2))&
!&          - (  ss(i,1)*Sr(1) - cs(i,1)*Si(1))&
!&       ) )
     __INCR( el_dz(i), efac * q(i) * kk * ( s(1) + s(3) - ( s(2) + s(4) ) ) )!&
!&           + ( -ss(i,1)*Sr(1) + cs(i,1)*Si(1))&
!&           - ( -ss(i,2)*Sr(2) + cs(i,2)*Si(2))&
!&           + ( -ss(i,3)*Sr(3) + cs(i,3)*Si(3))&
!&           - ( -ss(i,4)*Sr(4) + cs(i,4)*Si(4))&
!&          + (  ss(i,4)*Sr(4) - cs(i,4)*Si(4))&
!&          - (  ss(i,3)*Sr(3) - cs(i,3)*Si(3))&
!&          + (  ss(i,2)*Sr(2) - cs(i,2)*Si(2))&
!&          - (  ss(i,1)*Sr(1) - cs(i,1)*Si(1))&
!&       ) )
    enddo ! npt
   enddo !kn
  enddo !jn
 enddo !in
! normalize properly electrostatic potential and gradients ; no x2 here because el has 1/2 prefactor ;
 el=el * oeps * oV
 el_dx = 2d0 * el_dx * oeps * oV ! needs 2 because : (i) account for -kx wavenumbers (2) |S(k)|^2 derivative has a 2 prefactor ; 1/2 in front of el. pot.
 el_dy = 2d0 * el_dy * oeps * oV
 el_dz = 2d0 * el_dz * oeps * oV
!
 __FREE(ss)
 __FREE(cs)
end subroutine direct3
