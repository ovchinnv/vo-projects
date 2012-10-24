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
module rng
 implicit none
 logical :: random_initialized=.false.
 integer :: s(4)=(/1,2,3,4/)
 private s
 contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine random_init(seeds)
 use parser
 use output
!
 character(len=200) :: msg___(20)=(/'','','','','','','','','','','','','','','','','','','',''/); integer :: i_
 integer, optional :: seeds(4)
 if (random_initialized) return
!
 if (present(seeds)) then ; s=seeds ;
 else
!
  if (existtag_nocase('random_seeds')) then ;
   s=atofv(getval_nocase('random_seeds'),4);
  else
   write(msg___,*)'SEEDS NOT SPECIFIED, USING [',s,']';call warning('RANDOM_INIT', msg___(1), 0)
  endif
 endif ! seeds
!
 if (.not.fatal_warning()) then
  call clcginit(s)
  random_initialized=.true.
 else
  random_initialized=.false.
 endif
!
 end subroutine random_init
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine random_reinit(seeds)
 use parser
 use output
!
 character(len=200) :: msg___(20)=(/'','','','','','','','','','','','','','','','','','','',''/); integer :: i_
 integer, optional :: seeds(4)
!
 if (present(seeds)) then ; s=seeds ;
 else
!
  if (existtag_nocase('random_seeds')) then ;
   s=atofv(getval_nocase('random_seeds'),4);
  else
   write(msg___,*)'SEEDS NOT SPECIFIED, USING [',s,']';call warning('RANDOM_INIT', msg___(1), 0)
  endif
 endif ! seeds
!
 if (.not.fatal_warning()) then
  call clcginit(s)
  random_initialized=.true.
 endif
!
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
  call warning(whoami, 'RNG NOT INTIALIZED. ABORT',-1)
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
  call warning(whoami, 'RNG NOT INTIALIZED. ABORT',-1)
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
