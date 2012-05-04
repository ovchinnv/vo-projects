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
module rng
 implicit none
 logical :: random_initialized=.false.
 contains
 subroutine random_init(seeds)
 integer :: seeds(4)
 if (random_initialized) return
 call clcginit(seeds)
 random_initialized=.true.
 end subroutine random_init
!
 subroutine random_reinit(seeds)
 integer :: seeds
 call clcginit(seeds)
 random_initialized=.true.
 end subroutine random_reinit
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function randomu(channel)
!
! uniformly distributed random number
!
 use output
!
 real*8 :: randomu
 integer, optional :: channel
 integer :: chan
 character(len=8), parameter :: whoami='RANDOMU'
!
 real*8 :: random
!
 if (present(channel)) then ; chan=channel ; else ; chan=1 ; endif
 if (random_initialized) then
  randomu=random(chan)
 else
  call error(whoami, 'RNG NOT INTIALIZED. ABORT',-1)
  randomu = -1d0
 endif
!
 end function randomu
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine randomu_vector(a, n, channel)
!
! uniformly distributed random numbers
!
 use output
!
 real*8 :: a(n)
 integer :: n
 integer, optional :: channel
 integer :: chan
 character(len=14), parameter :: whoami='RANDOMU_VECTOR'
!
 real*8 :: random
 integer :: i
!
 if (present(channel)) then ; chan=channel ; else ; chan=1 ; endif
 if (random_initialized) then
!
  do i=1,n ; a(i)=random(chan) ; enddo
!
 else
  call error(whoami, 'RNG NOT INTIALIZED. ABORT',-1)
 endif
!
 end subroutine randomu_vector
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine randomg_vector(a, n, channel)
!
! generate Gaussian random numbers using Box-Muller scheme
!
 use constants
 use output
!
 real*8 :: a(n)
 integer :: n
 integer, optional :: channel
 integer :: chan
 character(len=14), parameter :: whoami='RANDOMG_VECTOR'
!
 real*8 :: random, u1, u2
 integer :: i, m
!
 if (present(channel)) then ; chan=channel ; else ; chan=1 ; endif
 if (random_initialized) then
!
  i=1;
  m=mod(n,2);
  do while (i.lt.(n-m))
     u1=sqrt(-2d0*log(random(chan)))
     u2=twopi*random(chan)
     a(i)=u1*cos(u2) ; i=i+1
     a(i)=u1*sin(u2) ; i=i+1
  enddo
  if (m.eq.1) then
     u1=sqrt(-2d0*log(random(chan)))
     u2=twopi*random(chan)
     a(i)=u1*cos(u2)
  endif
!
 else
  call error(whoami, 'RNG NOT INTIALIZED. ABORT',-1)
 endif
!
 end subroutine randomg_vector
end module rng
